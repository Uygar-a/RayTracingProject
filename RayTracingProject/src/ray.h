#pragma once

#include "../glm/glm.hpp"


class ray {
public:
    ray(const glm::vec4& origin, const glm::vec4& direction)
        : o(origin), d(direction)
    {}

    glm::vec4 origin() const { return o; }
    glm::vec4 direction() const { return d; }

    glm::vec4 at(float t) const {
        return o + t * d;
    }


public:
    glm::vec4 o;
    glm::vec4 d;
};
