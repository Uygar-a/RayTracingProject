#pragma once
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
//https://github.com/tinyobjloader/tinyobjloader

#include "../glm/glm.hpp"
#include "triangle.h"

class object3d {
public:
    object3d(std::string objname)
        : obj(objname)
    {}
    object3d(std::string objname ,std::vector<triangle> Triangles, std::vector<glm::vec4> TriDiffColor,
        std::vector<glm::vec4> TriSpecColor, glm::vec4 minPoint, glm::vec4 maxPoint)

        : obj(objname), tri(Triangles), diff(TriDiffColor), spec(TriSpecColor), min(minPoint), max(maxPoint)
    {}

    std::vector<triangle> Triangles() const {
        return tri;
    }
    std::vector<glm::vec4> TriDiffColor() const {
        return diff;
    }
    std::vector<glm::vec4> TriSpecColor() const {
        return spec; 
    }

    std::string objname() const {
        return obj;
    }

    glm::vec4 maxPoint() const {
        return max;
    }

    glm::vec4 minPoint() const {
        return min;
    }
    
    int mtlYesNo() const {
        return mtl_exist;
    }


    void SetTriangles(std::vector<triangle> Triangles){
        tri = Triangles;
    }
    void SetTriDiffColor(std::vector<glm::vec4> TriDiffColor) {
        diff = TriDiffColor;
    }

    void SetTriSpecColor(std::vector<glm::vec4> TriSpecColor) {
        spec = TriSpecColor;
    }

    void SetminPoint(glm::vec4 minPoint) {
        min = minPoint;
    }
    void SetmaxPoint(glm::vec4 maxPoint) {
        max = maxPoint;
    }

    void Setmtlexist(int istheremtl) {
        mtl_exist = istheremtl;
    }

public:
    glm::vec4 min, max;
    std::string obj;
    std::vector<triangle> tri;
    std::vector<glm::vec4> diff;
    std::vector<glm::vec4> spec;
    int mtl_exist;
};