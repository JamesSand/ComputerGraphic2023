#ifndef CAMERA_H
#define CAMERA_H

#include <float.h>

#include <cmath>

#include "ray.hpp"
#include "utils.hpp"
#include <vecmath.h>
const float INF_FOCAL_LENGTH = 0x3f3f3f3f;
class Camera {
   public:
    Camera(const Vector3f &center, const Vector3f &direction,
           const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    void setCenter(const Vector3f &pos) { this->center = pos; }
    Vector3f getCenter() const { return this->center; }

    void setRotation(const Matrix3f &mat) {
        this->horizontal = mat.getCol(0);
        this->up = -mat.getCol(1);
        this->direction = mat.getCol(2);
    }
    Matrix3f getRotation() const {
        return Matrix3f(this->horizontal, -this->up, this->direction);
    }

    virtual void resize(int w, int h) {
        width = w;
        height = h;
    }

   protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

class PerspectiveCamera : public Camera {
   public:
    float getFovy() const { return fovyd; }

    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, float angle,
                      float f = 20.0f, float aperture = 1.0f)
        : Camera(center, direction, up, imgW, imgH),
          focalLength(f),
          aperture(aperture) {
        // angle is fovy in radian.
        fovyd = angle / 3.1415 * 180.0;
        fx = fy = (float)height / (2 * tanf(angle / 2));
        cx = width / 2.0f;
        cy = height / 2.0f;
    }

    void resize(int w, int h) override {
        fx *= (float)h / height;
        fy = fx;
        Camera::resize(w, h);
        cx = width / 2.0f;
        cy = height / 2.0f;
    }

    Ray generateRay(const Vector2f &point) override {
        // 模拟光圈带来的扰动
        float symm_rand_x = symmatric_rand();
        float symm_rand_y = symmatric_rand();
        float dx = symm_rand_x * aperture;
        float dy = symm_rand_y * aperture;

        // 这里考虑到要对基于光圈扰动之后的的射线计算方向，所以要整个乘一个 focal
        // float csx = focalLength * (point.x() - cx) / fx;
        // float csy = focalLength * (point.y() - cy) / fy;
        // Vector3f noised_direction(csx - dx, -csy - dy, focalLength);

        float origin_x = (point.x() - cx) / fx;
        float origin_y = (cy - point.y()) / fy;
        Vector3f noised_direction(origin_x - dx / focalLength, origin_y - dy / focalLength, 1.0);
        // 向世界坐标系转换的矩阵 R
        Matrix3f R(horizontal, -up, direction);
        noised_direction = (R * noised_direction).normalized();
        // 这里考虑了景深，所以要根据光圈对 ray 的 origin做扰动，
        // 具体来说，是在成像点附近做了一个光圈的矩形，在这上边随机扰动
        Vector3f noised_center = center +  horizontal * dx - up * dy;
        Ray ray(noised_center, noised_direction);
        return ray;
    }

   protected:
    // Perspective intrinsics
    float fx;
    float fy;
    float cx;
    float cy;
    float fovyd;
    float aperture, focalLength;
};
#endif  // CAMERA_H
