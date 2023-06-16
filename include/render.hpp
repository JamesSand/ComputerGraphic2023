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
static Vector3f ptColor(Ray ray, const Scene& scene) {
    Group* group = scene.getGroup();
    int depth = 0;
    Vector3f color(0, 0, 0), cf(1, 1, 1);
    while (true) {
        if (++depth > TRACE_DEPTH || cf.max() < 1e-3) return color;
        // 判断camRay是否和场景有交点,返回最近交点的数据,存储在hit中.
        Hit hit;
        if (!group->intersect(ray, hit)) {
            color += scene.getBackgroundColor();
            return color;
        }

        // Path Tracing
        ray.origin += ray.direction * hit.t;
        Material* material = hit.material;
        Vector3f refColor(hit.color), N(hit.normal);

        // Emission
        color += material->emission * cf;
        cf = cf * refColor;
        float type = RND2;
        if (type <= material->type.x()) {  // diffuse
            ray.direction = diffDir(N);
        } else if (type <=
                   material->type.x() + material->type.y()) {  // specular
            float cost = Vector3f::dot(ray.direction, N);
            ray.direction = (ray.direction - N * (cost * 2)).normalized();
        } else {  // refraction
            float n = material->refr;
            float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
            if (Vector3f::dot(N, ray.direction) > 0) {  // inside the medium
                N.negate();
                n = 1 / n;
            }
            n = 1 / n;
            float cost1 = -Vector3f::dot(N, ray.direction);  // cosine theta_1
            float cost2 =
                1.0 - n * n * (1.0 - cost1 * cost1);  // cosine theta_2
            float Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1,
                                                5.0);  // Schlick-approximation
            if (cost2 > 0 && RND2 > Rprob) {           // refraction direction
                ray.direction =
                    ((ray.direction * n) + (N * (n * cost1 - sqrt(cost2))))
                        .normalized();
            } else {  // reflection direction
                ray.direction = (ray.direction + N * (cost1 * 2));
            }
        }
    }
}

static Vector3f rcColor(Ray ray, const Scene& scene) {
    Group* group = scene.getGroup();
    int depth = 0;
    Vector3f color(0, 0, 0);
    // 判断camRay是否和场景有交点,返回最近交点的数据,存储在hit中.
    Hit hit;
    if (!group->intersect(ray, hit)) {
        color += scene.getBackgroundColor();
        return color;
    }
    for (int li = 0; li < scene.getNumLights(); ++li) {
        Light* light = scene.getLight(li);
        Vector3f L, lightColor;
        // 获得光照强度
        light->getIllumination(ray.pointAtParameter(hit.getT()), L, lightColor);
        // 计算局部光强
        color += hit.getMaterial()->phongShade(ray, hit, L, lightColor);
    }
}

// class PathTracer {
//    public:
//     const Scene& scene;
//     int samps;
//     const char* fout;
//     Vector3f (*radiance)(Ray ray, const Scene& scene);
//     PathTracer(const Scene& scene, int samps, const char* method,
//                const char* fout)
//         : scene(scene), samps(samps), fout(fout) {
//         if (!strcmp(method, "rc"))
//             radiance = rcColor;
//         else if (!strcmp(method, "pt"))
//             radiance = ptColor;
//         else {
//             cout << "Unknown method: " << method << endl;
//             exit(1);
//         }
//     }

//     void render() {
//         Camera* camera = scene.getCamera();
//         int w = camera->getWidth(), h = camera->getHeight();
//         cout << "Width: " << w << " Height: " << h << endl;
//         Image outImg(w, h);
//         time_t start = time(NULL);
// #pragma omp parallel for schedule(dynamic, 1)  // OpenMP
//         for (int y = 0; y < h; ++y) {
//             float elapsed = (time(NULL) - start), progress = (1. + y) / h;
//             fprintf(stderr, "\rRendering (%d spp) %5.2f%% Time: %.2f/%.2f sec",
//                     samps, progress * 100., elapsed, elapsed / progress);
//             for (int x = 0; x < w; ++x) {
//                 Vector3f color = Vector3f::ZERO;
//                 for (int s = 0; s < samps; ++s) {
//                     Ray camRay =
//                         camera->generateRay(Vector2f(x + RND, y + RND));
//                     color += radiance(camRay, scene);
//                 }
//                 outImg.SetPixel(x, y, color / samps);
//             }
//         }
//         outImg.SaveBMP(fout);
//     }
// };

class SPPM {
   public:
    const Scene& scene;
    int numRounds, numPhotons, ckpt_interval;
    std::string outdir;
    int w, h;
    Camera* camera;
    vector<Hit*> hitPoints;
    HitKDTree* hitKDTree;
    vector<Object3D*> illuminants;
    Group* group;


    SPPM(const Scene& scene, int numRounds, int numPhotons, int ckpt,
         const char* dir)
        : scene(scene),
          numRounds(numRounds),
          numPhotons(numPhotons),
          ckpt_interval(ckpt),
          outdir(dir) {
        camera = scene.getCamera();
        group = scene.getGroup();
        illuminants = group->getIlluminant();

        w = camera->getWidth();
        h = camera->getHeight();
        hitKDTree = nullptr;
        for (int u = 0; u < w; ++u)
            for (int v = 0; v < h; ++v) hitPoints.push_back(new Hit());
        cout << "Width: " << w << " Height: " << h << endl;
    }

    ~SPPM() {
        for (int u = 0; u < w; ++u)
            for (int v = 0; v < h; ++v) delete hitPoints[u * w + v];
        delete hitKDTree;
    }

    void forward() {
#pragma omp parallel for schedule(dynamic, 1)

        for(int i = 0; i < w; i++){
            for(int j = 0; j < h; j++){
                
                Ray ray = camera->generateRay(Vector2f(i + RND, j + RND));
                int hit_index = i * h + j;
                Hit * hit = hitPoints[hit_index];

                int depth = 0;
                Vector3f attenuation(1, 1, 1);
                while (true) {
                    if (++depth > TRACE_DEPTH || attenuation.max() < 1e-3){
                        break;
                    }
                    hit->t = INF;
                    if (!group->intersect(ray, *hit)) {
                        hit->fluxLight += hit->attenuation*scene.getBackgroundColor();
                        break;
                    }
                    ray.origin += ray.direction * (*hit).t;
                    Material* material = (*hit).material;
                    Vector3f N(hit->normal);
                    float type = RND2;
                    if (type <= material->type.x()) {  // Diffuse
                        hit->attenuation = attenuation * hit->color;
                        hit->fluxLight += hit->attenuation * material->emission;
                        break;
                    } else if (type <= material->type.x() + material->type.y()) {
                        float cost = Vector3f::dot(ray.direction, N);
                        ray.direction = (ray.direction - N * (cost * 2)).normalized();
                    } else {
                        float n = material->refr;
                        float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
                        if (Vector3f::dot(N, ray.direction) > 0) {  // inside the medium
                            N.negate();
                            n = 1 / n;
                        }
                        n = 1 / n;
                        float cost1 =
                            -Vector3f::dot(N, ray.direction);  // cosine theta_1
                        float cost2 =
                            1.0 - n * n * (1.0 - cost1 * cost1);  // cosine theta_2
                        float Rprob =
                            R0 + (1.0 - R0) * pow(1.0 - cost1,
                                                5.0);   // Schlick-approximation
                        if (cost2 > 0 && RND2 > Rprob) {  // refraction direction
                            ray.direction =
                                ((ray.direction * n) + (N * (n * cost1 - sqrt(cost2))))
                                    .normalized();
                        } else {  // reflection direction
                            ray.direction =
                                (ray.direction + N * (cost1 * 2)).normalized();
                        }
                    }
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
                    Ray ray = illuminants[j]->randomRay(-1, (long long)round * numPhotons + (round + 1) * w * h + i);
                    
                    long long seed = round * numPhotons + i;

                    int depth = 0;
                    Vector3f attenuation = color * Vector3f(250, 250, 250);
                    while (true) {
                        if (++depth > TRACE_DEPTH || attenuation.max() < 1e-3){
                            break;
                        }
                        Hit hit;
                        if (!group->intersect(ray, hit)){
                            break;
                        }
                        ray.origin += ray.direction * hit.t;
                        Material* material = hit.material;
                        Vector3f N(hit.normal);
                        float type = RND2;
                        if (type <= material->type.x()) {  // Diffuse
                            hitKDTree->update(hitKDTree->root, hit.p, attenuation,
                                            ray.direction);
                            ray.direction = diffDir(N, -1, seed);
                        } else if (type <= material->type.x() + material->type.y()) {
                            float cost = Vector3f::dot(ray.direction, N);
                            ray.direction = (ray.direction - N * (cost * 2)).normalized();
                        } else {
                            float n = material->refr;
                            float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
                            if (Vector3f::dot(N, ray.direction) > 0) {  // inside the medium
                                N.negate();
                                n = 1 / n;
                            }
                            n = 1 / n;
                            float cost1 =
                                -Vector3f::dot(N, ray.direction);  // cosine theta_1
                            float cost2 =
                                1.0 - n * n * (1.0 - cost1 * cost1);  // cosine theta_2
                            float Rprob =
                                R0 + (1.0 - R0) * pow(1.0 - cost1,
                                                    5.0);   // Schlick-approximation
                            if (cost2 > 0 && RND2 > Rprob) {  // refraction direction
                                ray.direction =
                                    ((ray.direction * n) + (N * (n * cost1 - sqrt(cost2))))
                                        .normalized();
                            } else {  // reflection direction
                                ray.direction =
                                    (ray.direction + N * (cost1 * 2)).normalized();
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

            // 重建 kd tree
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