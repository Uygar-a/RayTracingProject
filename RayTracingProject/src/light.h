#pragma once

#include "../glm/glm.hpp"

#include <iostream>

class light {
public:
    light(const glm::vec4& origin, const glm::vec4& source_color)
        : o(origin), c(source_color)
    {}

    glm::vec4 origin() const { return o; }
    glm::vec4 lightcolor() const { return c; }

public:
    glm::vec4 o;
    glm::vec4 c;
};


