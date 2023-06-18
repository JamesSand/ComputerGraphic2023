#ifndef TEXTURE_H
#define TEXTURE_H
#include <string>

#include <vecmath.h>
using std::string;
class Texture {  // 纹理
    public:

    int counter;

    unsigned char *pic;
    int w, h, c;
    Texture(const char *textureFile);

    Vector3f getColor(float u, float v) const;
    Vector3f getColor(int idx) const;
    Vector3f getColor(int u, int v) const;
    float getDisturb(float u, float v, Vector2f &grad) const;
    inline int getIdx(float u, float v) const;
    inline float getGray(int idx) const;

    Vector3f get_color_mod(float u, float v) const;

};

#endif  // !TEXTURE_H