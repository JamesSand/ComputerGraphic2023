#ifndef OBJECTKDTREE_H
#define OBJECTKDTREE_H

#define MAX_TREE_DEPTH 24
#define MAX_FACE_NUM 128


#include <vecmath.h>

#include <map>
#include <vector>

#include "bound.hpp"
#include "hit.hpp"
#include "object3d.hpp"
using std::map;
using std::vector;
class ObjectKDTreeNode {
   public:
    Vector3f min, max;
    vector<Object3D*>* faces;
    ObjectKDTreeNode *ls, *rs;
    int l, r;
    bool inside(Object3D* face) {
        Vector3f faceMin = face->min();
        Vector3f faceMax = face->max();
        return (faceMin.x() < max.x() ||
                faceMin.x() == max.x() && faceMin.x() == faceMax.x()) &&
               (faceMax.x() > min.x() ||
                faceMax.x() == min.x() && faceMin.x() == faceMax.x()) &&
               (faceMin.y() < max.y() ||
                faceMin.y() == max.y() && faceMin.y() == faceMax.y()) &&
               (faceMax.y() > min.y() ||
                faceMax.y() == min.y() && faceMin.y() == faceMax.y()) &&
               (faceMin.z() < max.z() ||
                faceMin.z() == max.z() && faceMin.z() == faceMax.z()) &&
               (faceMax.z() > min.z() ||
                faceMax.z() == min.z() && faceMin.z() == faceMax.z());
    }
};

class ObjectKDTree {

    public:

    int n;
    Vector3f** vertices;

    ObjectKDTreeNode* root;
    vector<Object3D*>* faces;
    
    // //  ========================================= intersect part

    bool intersect(const Ray& ray, Hit& hit) const {
        Object3D* nextFace = nullptr;
        // start search from root
        return intersect(root, ray, nextFace, hit);
    }

    bool intersect(ObjectKDTreeNode* p, const Ray& ray, Object3D*& nextFace,
                   Hit& hit) const {
        bool flag = false;
        for (int i = 0; i < p->faces->size(); ++i)
            if ((*p->faces)[i]->intersect(ray, hit)) {
                nextFace = (*p->faces)[i];
                flag = true;
            }
        float tl = cuboidIntersect(p->ls, ray),
              tr = cuboidIntersect(p->rs, ray);
        if (tl < tr) {
            if (hit.t <= tl) return flag;
            if (p->ls) flag |= intersect(p->ls, ray, nextFace, hit);
            if (hit.t <= tr) return flag;
            if (p->rs) flag |= intersect(p->rs, ray, nextFace, hit);
        } else {
            if (hit.t <= tr) return flag;
            if (p->rs) flag |= intersect(p->rs, ray, nextFace, hit);
            if (hit.t <= tl) return flag;
            if (p->ls) flag |= intersect(p->ls, ray, nextFace, hit);
        }
        return flag;
    }


    float cuboidIntersect(ObjectKDTreeNode* p, const Ray& ray) const {
        float t = INF;
        if (!p) return t;
        AABB(p->min, p->max).intersect(ray, t);
        return t;
    }

    

    // //  ========================================= build tree part

    ObjectKDTree(vector<Object3D*>* input_faces) {
        Vector3f low = Vector3f(INF, INF, INF);
        Vector3f high = Vector3f(-INF, -INF, -INF);

        for (auto face : *input_faces) {
            low = min_Vf3(low, face->min());
            high = max_Vf3(high, face->max());
        }
        // start to build obj kd tree
        root = build_obj_kdtree(1, 0, input_faces, low, high);

        faces = new vector<Object3D*>;
        collect_face(root, faces);
    }

    ObjectKDTreeNode* build_obj_kdtree(int depth, int dimension, vector<Object3D*>* faces,
                        const Vector3f& low, const Vector3f& high) {
        ObjectKDTreeNode* current_node = new ObjectKDTreeNode;
        current_node->min = low;
        current_node->max = high;

        // 遍历父节点传进来的所有 face 
        current_node->faces = new vector<Object3D*>;
        for (auto face : *faces){
            if (current_node->inside(face)) {
                current_node->faces->push_back(face);
            }
        }

        if (current_node->faces->size() > MAX_FACE_NUM && depth < MAX_TREE_DEPTH) {

            // calculate next middle point accord to dimension
            float max_l_x = current_node->max.x();
            float max_l_y = current_node->max.y();
            float max_l_z = current_node->max.z();
            float min_r_x = current_node->min.x();
            float min_r_y = current_node->min.y();
            float min_r_z = current_node->min.z();

            if (dimension == 0){
                float mid_x = (current_node->min.x() + current_node->max.x()) / 2;
                max_l_x = mid_x;
                min_r_x = mid_x;
            }else if (dimension == 1){
                float mid_y = (current_node->min.y() + current_node->max.y()) / 2;
                max_l_y = mid_y;
                min_r_y = mid_y;
            }else{
                float mid_z = (current_node->min.z() + current_node->max.z()) / 2;
                max_l_z = mid_z;
                min_r_z = mid_z;
            }

            Vector3f maxL(max_l_x, max_l_y, max_l_z);
            Vector3f minR(min_r_x, min_r_y, min_r_z);

            // recursive build left and right
            current_node->ls = build_obj_kdtree(depth + 1, (dimension + 1) % 3, current_node->faces, low, maxL);
            current_node->rs = build_obj_kdtree(depth + 1, (dimension + 1) % 3, current_node->faces, minR, high);

            vector<Object3D*>*faceL = current_node->ls->faces, *faceR = current_node->rs->faces;
            
            // 在左右子树中都出现的面片放到父节点上
            map<Object3D*, int> face_counter;

            for (auto face : *faceL){
                face_counter[face]++;
            } 

            for (auto face : *faceR){
                face_counter[face]++;
            } 

            current_node->ls->faces = new vector<Object3D*>;
            current_node->rs->faces = new vector<Object3D*>;
            current_node->faces->clear();

            for (auto face : *faceL){
                if (face_counter[face] == 1){
                    current_node->ls->faces->push_back(face);
                }else{
                    current_node->faces->push_back(face);
                }
            }
                
            for (auto face : *faceR){
                if (face_counter[face] == 1){
                    current_node->rs->faces->push_back(face);
                }
            }
                
        } else{
            // 超过了最多面片数量或者递归深度，不再继续生长节点
            current_node->ls = nullptr;
            current_node->rs = nullptr;
        }
            
        return current_node;
    }

    void collect_face(ObjectKDTreeNode* current_node, vector<Object3D*>* faces) {
        current_node->l = faces->size();
        for (auto face : *(current_node->faces)){
            faces->push_back(face);
        } 
        current_node->r = faces->size();

        if (current_node->ls){
            collect_face(current_node->ls, faces);
        }
        if (current_node->rs){
            collect_face(current_node->rs, faces);
        } 
    }

    //  ========================================= utils function

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

};

#endif  // !OBJECTKDTREE_H