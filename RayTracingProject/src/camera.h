#pragma once

#include "ray.h"

//Camera library is implemented from https://raytracing.github.io/books/RayTracingInOneWeekend.html

class camera {
public:
    camera(
        glm::vec4 lookfrom,
        glm::vec4 lookat,
        glm::vec4 vup,
        float vfov, // vertical field-of-view in degrees
        float aspect_ratio
    ) {
        const double pi = 3.14159265359;

        float theta = vfov * pi/180.0;

        float h = tan(theta / 2);
        float viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        auto w = glm::normalize(lookfrom - lookat);
        auto u = glm::normalize(homCross(vup, w));
        auto v = homCross(w, u);

        origin = lookfrom;
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = origin - horizontal / 2.0f - vertical / 2.0f - w;
    }

    ray get_ray(float s, float t) const {
        return ray(origin, lower_left_corner + s * horizontal + t * vertical - origin);
    }

    glm::vec4 getOrigin() {
        return origin;
    }

    void setOrigin(glm::vec4 camOrigin) {
        origin = camOrigin;
    }


    //cross product for homogenous notation
    glm::vec4 homCross(glm::vec4 v1, glm::vec4 v2) {
        glm::vec3 vec1 = glm::vec3(v1.x, v1.y, v1.z);
        glm::vec3 vec2 = glm::vec3(v2.x, v2.y, v2.z);
        glm::vec3 crossed = glm::cross(vec1, vec2);

        return glm::vec4(crossed.x, crossed.y, crossed.z, 0);
    }

private:
    glm::vec4 origin;
    glm::vec4 lower_left_corner;
    glm::vec4 horizontal;
    glm::vec4 vertical;
};