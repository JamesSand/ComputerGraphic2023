#ifndef PATH_TRACER_H
#define PATH_TRACER_H

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "camera.hpp"
#include "constants.h"
#include "group.hpp"
#include "hit.hpp"
#include "hit_kdtree.hpp"
#include "image.hpp"
#include "light.hpp"
#include "ray.hpp"
#include "scene.hpp"
#include "utils.hpp"
using namespace std;

class SPPM {
   public:
    Scene * scene;
    int numRounds, numPhotons, ckpt_interval;
    std::string outdir;
    int w, h;
    Camera* camera;
    vector<Hit*> hitPoints;
    HitKDTree* hitKDTree;
    vector<Object3D*> illuminants;
    Group* group;


    SPPM(string input_file_path, int numRounds, int numPhotons, int ckpt,
         string dir)
        : 
          numRounds(numRounds),
          numPhotons(numPhotons),
          ckpt_interval(ckpt),
          outdir(dir) {

        scene = new Scene(input_file_path.c_str());
        camera = scene->getCamera();
        group = scene->getGroup();
        illuminants = group->getIlluminant();

        w = camera->getWidth();
        h = camera->getHeight();
        hitKDTree = nullptr;

        // 初始化 hitpoint
        for(int i = 0; i < w; i++){
            for(int j = 0; j < h; j++){
                hitPoints.push_back(new Hit());
            }
        }
        cout << "Image Width: " << w << " Height: " << h << endl;

    }

    ~SPPM() {
        for(int i = 0; i < w; i++){
            for(int j = 0; j < h; j++){
                int hit_index = i * h + j;
                delete hitPoints[hit_index];
            }
        }

        delete hitKDTree;
        delete scene;
    }

    void forward() {

#pragma omp parallel for schedule(dynamic, 1)

        for(int i = 0; i < w; i++){
            for(int j = 0; j < h; j++){
                
                // 抗锯齿
                float symm_rand_x = symmatric_rand();
                float symm_rand_y = symmatric_rand();
                
                Ray ray = camera->generateRay(Vector2f(i + symm_rand_x, j + symm_rand_y));
                int hit_index = i * h + j;
                Hit * hit = hitPoints[hit_index];

                int iter_counter = 0;
                // 初始光设为白光，之后的每次迭代，直接相乘即可
                Vector3f attenuation(1, 1, 1);
                while (true) {
                    iter_counter += 1;
                    if (iter_counter > MAX_ITERATION_NUM || attenuation.max() < STOP_ENERGY){
                        // 如果超过迭代次数，或者光线强度过低，则退出迭代
                        break;
                    }
                    // 现在是追踪camera 出来的光想，所以每次迭代需要重置 t
                    hit->t = INF;
                    if (!group->intersect(ray, *hit)) {
                        hit->fluxLight += hit->attenuation*scene->getBackgroundColor();
                        break;
                    }
                    // 每次迭代要将光线原点更新到交点
                    ray.origin += ray.direction * (*hit).t;
                    Material* material = (*hit).material;
                    Vector3f hit_normal(hit->normal);
                    // 根据随机数来决定是，漫反射，反射，还是折射
                    float rand_type = uniform_rand();
                    if (rand_type <= material->type.x()) {  
                        // 漫反射，吸收 break
                        // 吸收点颜色设置为当前光的颜色乘以反射面的颜色
                        hit->attenuation = attenuation * hit->color;

                        hit->fluxLight += hit->attenuation * material->emission;
                        break;
                        
                    } else if (rand_type <= material->type.x() + material->type.y()) {
                        // 镜面反射
                        // cos t 是入射方向与 hit normal 的夹角
                        float cost = Vector3f::dot(ray.direction, hit_normal);
                        ray.direction = (ray.direction - hit_normal * (cost * 2)).normalized();
                    } else {

                        // 如下代码的数学原理参考了
                        // https://zhuanlan.zhihu.com/p/375746359
                        // https://zhuanlan.zhihu.com/p/443186414

                        // 折射
                        float n = material->refr;
                        float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
                        if (Vector3f::dot(hit_normal, ray.direction) > 0) {  
                            // 如果从物体里向空气里折射，需要将 hit norm 反向，同时将折射率交换，也即去倒数
                            hit_normal.negate();
                            n = 1 / n;
                        }
                        n = 1 / n;
                        float cost1 = - Vector3f::dot(hit_normal, ray.direction);  
                        
                        // 这里计算 cost2 开根之前
                        // 考虑到有可能入射角太大，是镜面反射，不开方
                        float cost2_before_root = 1.0 - n * n * (1.0 - cost1 * cost1); 
                        
                        // 这一部分是菲涅尔方程的 Schlick 近似的结果
                        // 有 Rprob 的概率发生镜面反射
                        float Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); 
                        if (cost2_before_root > 0 && uniform_rand() > Rprob) { 
                            // 折射
                            ray.direction = ((ray.direction * n) + (hit_normal * (n * cost1 - sqrt(cost2_before_root)))).normalized();
                        } else {  // reflection direction
                            // 镜面反射
                            // 注意这里 hitnormal 的方向是负的
                            ray.direction = (ray.direction + hit_normal * (cost1 * 2)).normalized();
                        }
                    }
                    // 改变光的颜色
                    attenuation = attenuation * hit->color;
                }

            }
        }
    }

    void backward(int round) {
                int photonsPerLight = numPhotons / illuminants.size();

// photon tracing pass
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < photonsPerLight; ++i) {
                for (int j = 0;j < illuminants.size(); ++j) {
                
                    Vector3f color = illuminants[j]->material->emission;

                    long long rand_ray_seed = (long long)round * numPhotons + (round + 1) * w * h + i;

                    Ray ray = illuminants[j]->randomRay(-1, rand_ray_seed);
                    
                    long long seed = round * numPhotons + i;

                    int iter_count = 0;
                    // 这里或许是最后再乘比较好
                    Vector3f attenuation = color * Vector3f(250, 250, 250);
                    while (true) {
                        iter_count += 1;
                        if (iter_count > MAX_ITERATION_NUM || attenuation.max() < STOP_ENERGY){
                            break;
                        }

                        Hit hit;
                        if (!group->intersect(ray, hit)){
                            break;
                        }

                        ray.origin += ray.direction * hit.t;
                        Material* material = hit.material;
                        Vector3f hit_normal(hit.normal);
                        float type = uniform_rand();
                        if (type <= material->type.x()) {  // Diffuse
                            // 漫反射，要在 kdtree 里边记录
                            hitKDTree->update(hitKDTree->root, hit.p, attenuation,
                                            ray.direction);
                            ray.direction = diffDir(hit_normal, -1, seed);
                        } else if (type <= material->type.x() + material->type.y()) {
                            // 反射
                            float cost = Vector3f::dot(ray.direction, hit_normal);
                            ray.direction = (ray.direction - hit_normal * (cost * 2)).normalized();
                        } else {
                            // 折射
                            float n = material->refr;
                            float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
                            if (Vector3f::dot(hit_normal, ray.direction) > 0) {  // inside the medium
                                hit_normal.negate();
                                n = 1 / n;
                            }
                            n = 1 / n;
                            float cost1 = -Vector3f::dot(hit_normal, ray.direction);  // cosine theta_1
                            float cost2_before_root = 1.0 - n * n * (1.0 - cost1 * cost1);  // cosine theta_2
                            
                            float Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0);   // Schlick-approximation
                            if (cost2_before_root > 0 && uniform_rand() > Rprob) {  
                                // 折射
                                ray.direction = ((ray.direction * n) + (hit_normal * (n * cost1 - sqrt(cost2_before_root)))).normalized();
                            } else {  
                                // 反射
                                ray.direction = (ray.direction + hit_normal * (cost1 * 2)).normalized();
                            }
                        }
                        attenuation = attenuation * hit.color;
                    }

                } 
            }
    }

    void render() {
        // time_t start = time(NULL);
        time_t start_time = time(NULL);
        // Vector3f color = Vector3f::ZERO;
        for (int round = 0; round < numRounds; round++) {

            float time_cost = (time(NULL) - start_time);
            float done_percentage = float(round + 1) / numRounds;
            cout << "current iter " << round + 1 << " total iter " << numRounds;
            cout << " time elapse " << time_cost << " sec ";
            cout << "finish time " << time_cost / (round + 1) * numRounds << "sec";
            cout << endl;

            // 相机采集光子 forward
            forward();

            // 根据从 camera 发出的光线，重建 kd tree
            if (hitKDTree){
                delete hitKDTree;
            } 
            hitKDTree = new HitKDTree(&hitPoints);

            // 光子散落 backward
            backward(round);

            if ((round + 1) % ckpt_interval == 0) {
                Image ckpt_img(w, h);

                for(int i = 0; i < w; i++){
                    for (int j = 0; j < h; j++){
                        int hitpoint_index = i * h + j;
                        Hit * current_hit = hitPoints[hitpoint_index];

                        Vector3f current_color = current_hit->flux / (M_PI * current_hit->r2 * numPhotons * (round + 1)) +
                                                    current_hit->fluxLight / (round + 1);
                        ckpt_img.SetPixel(i, j, current_color);
                    }
                }

                // save checkpoint image
                char filename[100];
                sprintf(filename, "ckpt-%d.bmp", round + 1);
                // 这个地方很奇怪，似乎只能这么写才能正确写入
                ckpt_img.SaveBMP((outdir + "/" + filename).c_str());
                cout << "ckpt at iter " << round + 1 << endl;

            }
        }


        Image final_img(w, h);
        for(int i = 0; i < w; i++){
            for (int j = 0; j < h; j++){
                int hitpoint_index = i * h + j;
                Hit * current_hit = hitPoints[hitpoint_index];
                Vector3f current_color = current_hit->flux / (M_PI * current_hit->r2 * numPhotons * numRounds) +
                                            current_hit->fluxLight / numRounds;
                final_img.SetPixel(i, j, current_color);
            }
        }
        // save checkpoint image
        final_img.SaveBMP((outdir + "/" + "final_image.bmp").c_str());
        cout << "final image saved" << endl;

    }

};

#endif  // !PATH_TRACER_H