#ifndef HIT_KDTREE_H
#define HIT_KDTREE_H
#include <algorithm>

#include "constants.h"
#include "hit.hpp"

struct HitKDTreeNode {
//    public:
    Hit *hit;
    Vector3f min, max;
    float maxr2;
    HitKDTreeNode *ls, *rs;
};

// 这里只是沿用了 object kdtree 的写法，将 hitpoint 按照三维进行分类
class HitKDTree {

    public:

    int total_hits;
    Hit **hits;

    HitKDTreeNode *root;


    // =================== build_kdtree tree =====================


    HitKDTree(vector<Hit *> *input_hits) {
        total_hits = input_hits->size();

        hits = new Hit *[total_hits];
        for (int i = 0; i < total_hits; i++){
            hits[i] = (*input_hits)[i];
        }

        root = build_kdtree(0, total_hits - 1, 0);
    }

    ~HitKDTree() {
        if (root == nullptr){
            return;
        }

        delete_node(root);
        delete[] hits;
    }


    HitKDTreeNode *build_kdtree(int left, int right, int dimension) {

        HitKDTreeNode *current_node = new HitKDTreeNode;
        Vector3f temp_min = Vector3f(INF, INF, INF);
        Vector3f temp_max = Vector3f(-INF, -INF, -INF);
        float temp_radius = 0.0;

        for (int i = left; i <= right; i++) {
            temp_min = min_Vf3(temp_min, hits[i]->position);
            temp_max = max_Vf3(temp_max, hits[i]->position);
            temp_radius = max_float(temp_radius, hits[i]->r2);
        }

        current_node->min = temp_min;
        current_node->max = temp_max;
        current_node->maxr2 = temp_radius;

        int middle = (left + right) / 2;
        
        // set the middle to the right place
        if (dimension == 0){
            std::nth_element(hits + left, hits + middle, hits + right + 1, x_compare);
        }else if (dimension == 1){
            std::nth_element(hits + left, hits + middle, hits + right + 1, y_compare);
        }else{
            std::nth_element(hits + left, hits + middle, hits + right + 1, z_compare);
        }
            
        current_node->hit = hits[middle];

        if (left <= middle - 1){
            // build left
            current_node->ls = build_kdtree(left, middle - 1, (dimension + 1) % 3);
        }else{
            current_node->ls = nullptr;
        }
        

        if (middle + 1 <= right){
            // build right
            current_node->rs = build_kdtree(middle + 1, right, (dimension + 1) % 3);
        }else{
            current_node->rs = nullptr;
        }
            
        
        return current_node;
    }


    // =================== delete and update =====================


    bool aabb_acc_check(HitKDTreeNode *current_node, const Vector3f& position){
        float upper_bound = current_node->maxr2;
        float radius_counter = 0;

        Vector3f temp_max = current_node->max;
        Vector3f temp_min = current_node->min;

        if (position.x() > temp_max.x()) {
            radius_counter += square_float(position.x() - temp_max.x());
        }else if (position.x() < temp_min.x()){
            radius_counter += square_float(temp_min.x() - position.x());
        } 

        if (position.y() > temp_max.y()) {
            radius_counter += square_float(position.y() - temp_max.y());
        }else if (position.y() < temp_min.y()){
            radius_counter += square_float(temp_min.y() - position.y());
        } 

        if (position.z() > temp_max.z()) {
            radius_counter += square_float(position.z() - temp_max.z());
        }else if (position.z() < temp_min.z()){
            radius_counter += square_float(temp_min.z() - position.z());
        } 

        
        if (radius_counter > current_node->maxr2){
            return false;
        }

        return true;


    }

    void update(HitKDTreeNode *current_node, const Vector3f &photon,
                const Vector3f &attenuation, const Vector3f &direction) {

        // recursive base
        if (current_node == nullptr){
            return;
        }

        // aabb 包围盒加速检查
        if (aabb_acc_check(current_node, photon) == false){
            return;
        }

        if ((photon - current_node->hit->position).squaredLength() <= current_node->hit->r2) {
            // 小于 吸收 radius 可以被吸收
            Hit * hitpoint = current_node->hit;
            // 和之前做一个加权平均
            float factor = (hitpoint->phono_num * ALPHA + ALPHA) / (hitpoint->phono_num * ALPHA + 1.);
            // Vector3f dr = direction - hitpoint->normal * (2 * Vector3f::dot(direction, hitpoint->normal));
            hitpoint->phono_num++;
            // 本来 alpha 是 0.6 但是这里更新的是 radius square
            // 更新半径
            hitpoint->r2 *= factor;
            // 更新光通量
            hitpoint->flux = (hitpoint->flux + hitpoint->attenuation * attenuation) * factor;
        }

        // recursive update
        if (current_node->ls) {
            update(current_node->ls, photon, attenuation, direction);
        }
        if (current_node->rs){
            update(current_node->rs, photon, attenuation, direction);
        }

        // update max radius squre
        current_node->maxr2 = current_node->hit->r2;

        if (current_node->ls && current_node->ls->hit->r2 > current_node->maxr2){
            current_node->maxr2 = current_node->ls->hit->r2;
        }
        if (current_node->rs && current_node->rs->hit->r2 > current_node->maxr2){
            current_node->maxr2 = current_node->rs->hit->r2;
        }
    }
    
    void delete_node(HitKDTreeNode *current_node) {
        // recursive delete node
        if (current_node->ls){
            delete_node(current_node->ls);
        }
        if (current_node->rs){
            delete_node(current_node->rs);
        }
        delete current_node;
    }

   
    
    // ================================ utils ================================


    inline float square_float(float a) { return a * a; }

    inline float min_float(float a, float b){
        return a < b ? a : b;
    }

    inline float max_float(float a, float b){
        return a < b ? b : a;
    }

    inline Vector3f min_Vf3(const Vector3f& a, const Vector3f&b){
        return Vector3f(min_float(a.x(), b.x()),
                        min_float(a.y(), b.y()),
                        min_float(a.z(), b.z())
                        );
    }

    inline Vector3f max_Vf3(const Vector3f& a, const Vector3f&b){
        return Vector3f(max_float(a.x(), b.x()),
                        max_float(a.y(), b.y()),
                        max_float(a.z(), b.z())
                        );
    }


    static bool x_compare(Hit *a, Hit *b) { 
        return a->position.x() < b->position.x(); 
    }
    static bool y_compare(Hit *a, Hit *b) { 
        return a->position.y() < b->position.y(); 
    }
    static bool z_compare(Hit *a, Hit *b) { 
        return a->position.z() < b->position.z(); 
    }

};
#endif