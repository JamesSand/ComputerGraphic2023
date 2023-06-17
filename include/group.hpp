#ifndef GROUP_H
#define GROUP_H

#include <iostream>
#include <vector>

#include "hit.hpp"
#include "object3d.hpp"
#include "object_kdtree.hpp"
#include "ray.hpp"

class Group {
   public:
    Group(const vector<Object3D *> &objs) {
        for (auto obj : objs) {
            vector<Object3D *> temp_face = obj->getFaces();
            faces.insert(faces.end(), temp_face.begin(), temp_face.end());
        }
        obj_kdtree = new ObjectKDTree(&faces);
    }
    ~Group() { 
        delete obj_kdtree; 
    }

    // 对列表里所有物体都求一遍交点
    bool intersect(const Ray &r, Hit &h) { 
        return obj_kdtree->intersect(r, h); 
    }

    // bool sequentialSearch(const Ray &r, Hit &h) {
    //     bool flag = false;
    //     for (auto face : faces)
    //         if (face) flag |= face->intersect(r, h);
    //     return flag;
    // }

    int getGroupSize() { return faces.size(); }
    Object3D *operator[](const int &i) {
        if (i >= faces.size() || i < 0) {
            std::cout << "Invalid index " << i << std::endl;
            return nullptr;
        }

        return faces[i];
    }

    vector<Object3D *> getIlluminant() const {
        vector<Object3D *> illuminants;
        for (int i = 0; i < faces.size(); i++){
            if (faces[i]->material->emission != Vector3f::ZERO){
                illuminants.push_back(faces[i]);
            }
        }
        return illuminants;
    }

   private:
    ObjectKDTree *obj_kdtree;
    vector<Object3D *> faces;
};

#endif
