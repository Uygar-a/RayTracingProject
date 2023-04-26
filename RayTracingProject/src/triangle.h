#pragma once
#include "../glm/glm.hpp"

class triangle{
public:
    triangle(const glm::vec4& p1, const glm::vec4& p2, const glm::vec4& p3)
        : point1(p1), point2(p2), point3(p3)
    {}

    glm::vec4 p1() const { return point1;
    }
    glm::vec4 p2() const { return point2;
    }
    glm::vec4 p3() const { return point3; }


public:
    glm::vec4 point1, point2, point3;
    glm::vec4 e[3];
};