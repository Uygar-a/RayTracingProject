#include "ray.h"
#include "light.h"
#include "triangle.h"
#include "../glm/glm.hpp"
#include "../glm/gtx/transform.hpp"

#include "object3d.h"
#include "camera.h"

#include <time.h>
#include <iostream>
#include <fstream>  
#include <vector>
#include <sstream>

//cross product for homogenous notation
glm::vec4 homCross(glm::vec4 v1, glm::vec4 v2) {
    glm::vec3 vec1 = glm::vec3(v1.x, v1.y, v1.z);
    glm::vec3 vec2 = glm::vec3(v2.x, v2.y, v2.z);
    glm::vec3 crossed = glm::cross(vec1, vec2);

    return glm::vec4(crossed.x, crossed.y, crossed.z, 0);
}

//Axis Aligned Bounding Box
bool aabb(ray& r, glm::vec4& left_bottom, glm::vec4& top_right) {
    double tmin = 0;
    double tmax = 0;
    double temp;

    double txmin = (left_bottom.x - r.origin().x) / r.direction().x;
    double txmax = (top_right.x - r.origin().x) / r.direction().x;

    if (txmin > txmax) {
        temp = txmin;
        txmin = txmax;
        txmax = temp;
    }
    double tymin = (left_bottom.y - r.origin().y) / r.direction().y;
    double tymax = (top_right.y - r.origin().y) / r.direction().y;

    if (tymin > tymax) {
        temp = tymin;
        tymin = tymax;
        tymax = temp;
    }
    double tzmin = (left_bottom.z - r.origin().z) / r.direction().z;
    double tzmax = (top_right.z - r.origin().z) / r.direction().z;


    if (tzmin > tzmax) {
        temp = tzmin;
        tzmin = tzmax;
        tzmax = temp;
    }
    //std::cout << ty1 << std::endl;

    if (tymax > txmax) {
        tmax = txmax;
    }
    else {
        tmax = tymax;
    }

    if (tmax > tzmax) {
        tmax = tzmax;
    }

    if (tymin < txmin) {
        tmin = txmin;
    }
    else {
        tmin = tymin;
    }

    if (tmin < tzmin) {
        tmin = tzmin;
    }

    //if(tmax >= tmin)
    //std::cout << "TRUE" << std::endl;

    return tmax >= tmin;
}

//The intersection test I wrote from the slides in course
bool intersectRayTriangleBarycentric(const glm::vec4& orig, const glm::vec4& d,
    const glm::vec4& p0, const glm::vec4& p1, const glm::vec4& p2,
    float& t)
{
    glm::vec4 e1 = p2 - p0;
    glm::vec4 e2 = p1 - p0;
    glm::vec4 s = orig - p0;


    // check if parallel ?
    glm::vec4 crosslength = homCross(d, e2);

    //std::cout << crosslength.x <<" "<< crosslength.y << " " << crosslength.z << " " << std::endl;

    if (glm::length(crosslength) == 0)  //almost 0 
        return false;  //they are parallel so they don't intersect ! 

    // check if no edge
    if (glm::length(e1) == 0)  //almost 0 
        return false;

    float divider = 1 / glm::dot(homCross(d, e2), e1);

    t = divider * glm::dot(homCross(s, e1), e2);
    if (t < 0) return false;
    float b1 = divider * glm::dot(homCross(d, e2), s);
    float b2 = divider * glm::dot(homCross(s, e1), d);



    if (b1 >= 0 && b2 >= 0 && b1 + b2 <= 1 && t > 0)
        return true;
    else
        return false;
}



//Phong Model only for triangle surfaces
glm::vec4 phongIlluminationonTriangle(const ray& viewer, const light& pointlight, const glm::vec4& v0, const glm::vec4& v1, const glm::vec4& v2,
    const glm::vec4 surfacecolor, glm::vec4 ks, const glm::vec4& surfacepoint) {

    float pi = 3.141593;

    glm::vec4  e1 = v1 - v0;
    glm::vec4 e2 = v2 - v0;

    glm::vec4 N = homCross(e1, e2);  //N 
    N = glm::normalize(N);
    glm::vec4 l = pointlight.origin() - surfacepoint;
    l = glm::normalize(l);


    auto dot_product = glm::dot(N, l);

    if (dot_product < 0)
        //std::cerr << "N*L is NEGATIVE " << std::endl;
        dot_product = -dot_product;

    glm::vec4 diffuse_reflection = surfacecolor * (pointlight.lightcolor() * dot_product);
    //std::cout << surfacecolor << "  " << std::endl;

    glm::vec4 diffuse_reflection_tocamera = 3 * (1 / pi) * diffuse_reflection;
    //std::cout << diffuse_reflection_tocamera << "  " << std::endl;

    //specularreflection starts here

    glm::vec4 r = 2.0f * (l * N) * N - l;
    glm::vec4 v = viewer.origin() - surfacepoint;
    v = glm::normalize(v);
    r = glm::normalize(r);

    glm::vec4 specular_reflection = glm::vec4(1, 1, 1, 1) * pointlight.lightcolor() * dot_product
        * static_cast<float>(pow(fmax(dot(r, v), 0), 1));

    glm::vec4 specular_reflection_wthks = ks * specular_reflection;


    // std::cout << diffuse_reflection_tocamera << "  " << std::endl;

    glm::vec4 total_reflection = diffuse_reflection_tocamera + specular_reflection_wthks;

    if (total_reflection.x > 1)
        total_reflection.x = 1;
    if (total_reflection.y > 1)
        total_reflection.y = 1;
    if (total_reflection.z > 1)
        total_reflection.z = 1;

    return  total_reflection;
}


//Check to see if there is any object covering the path from lightsource to camerascene intersection
glm::vec4 shadowTest(const ray& r, light sourcelight, std::vector<triangle> Triangles, glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, float t_ray, glm::vec4 pixel_color) {
    float light_t;
    //std::cout << Triangles.size() << std::endl;
    for (int newi = 0; newi < Triangles.size(); newi++) {

        //If it is the same triangle
        if (v0.x == Triangles[newi].p1().x && v0.y == Triangles[newi].p1().y && v0.z == Triangles[newi].p1().z &&
            v1.x == Triangles[newi].p2().x && v1.y == Triangles[newi].p2().y && v1.z == Triangles[newi].p2().z &&
            v2.x == Triangles[newi].p3().x && v2.y == Triangles[newi].p3().y && v2.z == Triangles[newi].p3().z)
        {
        }
        else if (intersectRayTriangleBarycentric(r.at(t_ray), glm::normalize(sourcelight.origin() - r.at(t_ray)), Triangles[newi].p1(),
            Triangles[newi].p2(), Triangles[newi].p3(), light_t)) {
            return glm::vec4(0, 0, 0, 1);
        }
    }
    return pixel_color;
}


inline bool exists_test(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

int main() {

    clock_t tStart = clock();

    //Load an OBJ

    //std::vector<object3d> objects;

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string error;
    std::string warn;

    int mtlfile_exist;

    //AABB
    std::vector<glm::vec4> vertlist;
    glm::vec4 vertices;
    glm::vec4 largest;
    glm::vec4 smallest;

    float shift_by_x = 0;
    float shift_by_y = 0;


    //triagle
    std::vector<std::vector<triangle>> allTriangles;
    std::vector<std::vector<glm::vec4>> allTriSpecColor;
    std::vector<std::vector<glm::vec4>> allTriDiffColor;


    std::vector<triangle> Triangles;
    std::vector<glm::vec4> TriDiffColor;
    std::vector<glm::vec4> TriSpecColor;

    std::vector<object3d> allobjects;



    double distancing_constant = 5; //move away the objects from the camera
    auto focal_length = 1.0;

    //Hold multiple objects


    int maxTrianglessize = -1;

    std::string objname = "cube";
    std::string objname2 = "poly_dachshund";
    std::string objname3 = "gray+sphere";
    std::string objname4 = "sphere";


    object3d cube1(objname4);
    //object3d cube2(objname2);
    //object3d cube3(objname3);

    //object3d cube4(objname4);


    allobjects.emplace_back(cube1);
    //allobjects.emplace_back(cube2);
    //allobjects.emplace_back(cube3);
    //allobjects.emplace_back(cube4);


    for (int o = 0; o < allobjects.size(); o++)
    {
        std::string filename = "models/" + allobjects[o].objname() + ".obj";

        bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &error, filename.c_str(), "models");

        std::cout << materials.size();

        if (!warn.empty()) {
            std::cout << "WARN: " << warn << std::endl;
        }

        mtlfile_exist = 0;
        allobjects[o].Setmtlexist(0);

        if (exists_test("models/" + allobjects[o].objname() + ".mtl")) {
            mtlfile_exist = 1;
            allobjects[o].Setmtlexist(1);
        }

        if (!error.empty()) {
            std::cerr << "ERR: " << error << std::endl;
        }

        if (!ret) {
            printf("Failed to load/parse .obj.\n");
            return false;
        }

        double size_multiplier = 1;
        size_multiplier = size_multiplier / 1;


        // Adding vertices to a vector and finding largest and smallest values in vertices to build aabb
        double largest_x = -INFINITY;
        double largest_y = -INFINITY;
        double largest_z = -INFINITY;
        double smallest_x = INFINITY;
        double smallest_y = INFINITY;
        double smallest_z = INFINITY;

        for (size_t vec_start = 0; vec_start < attrib.vertices.size(); vec_start += 3) {

            glm::vec4 vertices(attrib.vertices[vec_start] * size_multiplier, attrib.vertices[vec_start + 1] * size_multiplier, attrib.vertices[vec_start + 2] * size_multiplier, 1);

            //FOR BOUNDING THE OBJECT
            if (vertices.x > largest_x)
                largest_x = vertices.x;
            if (vertices.y > largest_y)
                largest_y = vertices.y;
            if (vertices.z > largest_z)
                largest_z = vertices.z;

            if (vertices.x < smallest_x)
                smallest_x = vertices.x;
            if (vertices.y < smallest_y)
                smallest_y = vertices.y;
            if (vertices.z < smallest_z)
                smallest_z = vertices.z;

            //Resize the objects to have the similar sizes
            if (vec_start >= (attrib.vertices.size() - 4) && (largest_y - smallest_y > 2.2 || largest_y - smallest_y < 1.8)) {

                std::cout << largest_y - smallest_y << std::endl;

                size_multiplier = 2.0 / (largest_y - smallest_y);

                vec_start = 0;
                vertlist.clear();

                largest_x = -INFINITY;
                largest_y = -INFINITY;
                largest_z = -INFINITY;
                smallest_x = INFINITY;
                smallest_y = INFINITY;
                smallest_z = INFINITY;
            }

            vertlist.emplace_back(vertices);

            //std::cout << vertlist.at(somewherelike1or2) << std::endl;
        }


        largest = glm::vec4(largest_x, largest_y, largest_z, 1);
        smallest = glm::vec4(smallest_x, smallest_y, smallest_z, 1);
        std::cout << smallest.x << " " << smallest.y << " " << smallest.z << "SMALLEST" << std::endl;
        std::cout << largest.x << " " << largest.y << " " << largest.z << "LARGEST" << std::endl;

        if (o > 0)
            shift_by_x = largest_x - allobjects[o - 1].minPoint().x + 1.5;

        std::cout << shift_by_x << std::endl;

        shift_by_y = smallest_y;

        largest = largest - glm::vec4(shift_by_x, shift_by_y, distancing_constant * focal_length, 1);
        smallest = smallest - glm::vec4(shift_by_x, shift_by_y, distancing_constant * focal_length, 1);

        //Put values from tinyobjloader to a simplified vectors

        for (size_t i = 0; i < shapes.size(); i++) {

            size_t index_offset = 0;

            // For each face
            for (size_t f = 0; f < shapes[i].mesh.num_face_vertices.size(); f++) {
                size_t fnum = shapes[i].mesh.num_face_vertices[f];

                // For each vertex in the face
                // 
                tinyobj::index_t idx0 = shapes[i].mesh.indices[index_offset + 0];
                tinyobj::index_t idx1 = shapes[i].mesh.indices[index_offset + 1];
                tinyobj::index_t idx2 = shapes[i].mesh.indices[index_offset + 2]; // We may need to add normal and texcoord later.
                idx0.vertex_index;
                idx1.vertex_index;
                idx2.vertex_index;

                //Add 3 vertices of triangles point1, point2, point3
                Triangles.emplace_back(vertlist.at(idx0.vertex_index) - glm::vec4(shift_by_x, shift_by_y, distancing_constant * focal_length, 1),
                    vertlist.at(idx1.vertex_index) - glm::vec4(shift_by_x, shift_by_y, distancing_constant * focal_length, 1),
                    vertlist.at(idx2.vertex_index) - glm::vec4(shift_by_x, shift_by_y, distancing_constant * focal_length, 1));

                //Use mtl file if it exists
                if (mtlfile_exist == 1) {
                    TriDiffColor.emplace_back(glm::vec4(materials[shapes[i].mesh.material_ids[f]].diffuse[0],
                        materials[shapes[i].mesh.material_ids[f]].diffuse[1],
                        materials[shapes[i].mesh.material_ids[f]].diffuse[2], 1));

                    TriSpecColor.emplace_back(glm::vec4(materials[shapes[i].mesh.material_ids[f]].specular[0],
                        materials[shapes[i].mesh.material_ids[f]].specular[1],
                        materials[shapes[i].mesh.material_ids[f]].specular[2], 1));
                }
                index_offset += fnum;
            }
        }

        allTriangles.emplace_back(Triangles);
        allTriDiffColor.emplace_back(TriDiffColor);
        allTriSpecColor.emplace_back(TriSpecColor);

        allobjects[o].SetmaxPoint(largest);
        allobjects[o].SetminPoint(smallest);


        shapes.clear();
        materials.clear();
        vertlist.clear();

        Triangles.clear();
        TriDiffColor.clear();
        TriSpecColor.clear();

    }
    // Image
    const auto aspect_ratio = 4.0 / 3.0;
    const int image_width = 480;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Ccamera
    auto origin = glm::vec4(0, 2, 0, 1);

    camera cam(origin, glm::vec4(0, 2, -distancing_constant, 1), glm::vec4(0, 1, 0, 1), 90, aspect_ratio);

    //ray direction towards lower_left_corner 
    //auto lower_left_corner = origin - horizontal / 2.0f - vertical / 2.0f - glm::vec4(0, 0, focal_length, 0);  //-1.77778 -1 -1


    glm::vec4 v0(0, 2, 2, 1), v1(0, 0, 2, 1), v2(2, 2, 2, 1);
    v0 = v0 - glm::vec4(0, 0, 5 * focal_length, 1);
    v1 = v1 - glm::vec4(0, 0, 5 * focal_length, 1);
    v2 = v2 - glm::vec4(0, 0, 5 * focal_length, 1);

    light sourcelight(glm::vec4(10, 10, -2, 1), glm::vec4(1, 1, 1, 1));


    float t = 0;;

    // Render
    //Sending only rays to create an empty simple image is implemented from
    //https://raytracing.github.io/books/RayTracingInOneWeekend.html

    //Rotation part starts here
    int maximgno = 1;
    float pi = 3.14159265359;
    float rotationdegree = 360 / maximgno;
    std::string path = "sequence/";

    for (int imgno = 0; imgno < maximgno; imgno++) {

        float rotationradians = 2 * imgno * pi * (rotationdegree / 360);

        glm::mat4 translateMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, distancing_constant));
        //printf(" %f,  %f,  %f\n", (translateMatrix* cam.getOrigin()).x, (translateMatrix* cam.getOrigin()).y, (translateMatrix* cam.getOrigin()).z);

        glm::vec4 translated = translateMatrix * cam.getOrigin();

        glm::mat4 rotateMatrix = glm::rotate(glm::mat4(1.0f), rotationradians, glm::vec3(0.0f, 1.0f, 0.0f));

        glm::vec4 rotated = rotateMatrix * translated;
        glm::mat4 translateMatrix2 = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -distancing_constant));

        glm::vec4 translate2back = translateMatrix2 * rotated;

        camera cam(translate2back, glm::vec4(0, 2, -distancing_constant, 1), glm::vec4(0, 1, 0, 1), 90, aspect_ratio);

        //Rotation part ends here

        std::string file = path + std::to_string(imgno) + ".ppm";
        //std::ofstream outfile(file);
        std::ofstream outfile("output.ppm");
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";


        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {

                auto u = float(i) / (image_width - 1);
                auto v = float(j) / (image_height - 1);
                //send rays from lower left to top right 

                ray r = cam.get_ray(u, v);

                glm::vec4 pixel_color;
                float t_ray = INFINITY;
                float oldt_ray = t_ray;
                glm::vec4 withphong;

                //Defining default color in case mtl file does not exist
                glm::vec4 diffuse_color = glm::vec4(1, 0, 1, 1);
                glm::vec4 specular_color = glm::vec4(0, 0, 0, 1);


                pixel_color = glm::vec4(0.1, 0.1, 0.1, 1);


                //Add 0-level ground 
                glm::vec4 tg0(-8, 0, -2, 1);
                glm::vec4 tg1(8, 0, -12, 1);
                glm::vec4 tg2(8, 0, -2, 1);
                glm::vec4 tg3(-8, 0, -12, 1);
                glm::vec4 ground_color = glm::vec4(0, 1, 0, 1);
                glm::vec4 g0, g1, g2;

                //Two Triangles for ground

                for (int gr = 0; gr < 2; gr++) {
                    if (gr == 0) {
                        g0 = tg0; g1 = tg1; g2 = tg2;
                    }
                    else {
                        g0 = tg0; g1 = tg1; g2 = tg3;
                    }

                    if (intersectRayTriangleBarycentric(r.origin(), r.direction(), g0, g1, g2, t_ray)) {
                        if (t_ray < oldt_ray) {
                            ray shadow(r.at(t_ray), glm::normalize(sourcelight.origin() - r.at(t_ray)));

                            oldt_ray = t_ray;

                            pixel_color = phongIlluminationonTriangle(r, sourcelight, g0, g1, g2,
                                ground_color, specular_color, r.at(t_ray));

                            for (int o = 0; o < allobjects.size(); o++) {
                                //AA Bounding BOX parameters
                                auto left = allobjects[o].minPoint();
                                auto right = allobjects[o].maxPoint();
                                if (aabb(shadow, left, right))
                                    //Shadow on ground
                                    pixel_color = shadowTest(r, sourcelight, allTriangles[o], g0, g1, g2, t_ray, pixel_color);
                            }
                        }
                    }
                }

                for (int o = 0; o < allobjects.size(); o++) {

                    //AA Bounding BOX parameters
                    auto left = allobjects[o].minPoint();
                    auto right = allobjects[o].maxPoint();
                    //AA Bounding BOX
                    if (aabb(r, left, right)) {
                        //std::cout << right.y<<std::endl;

                        //take the vertices from vertlist using indices vector
                        //Triangles.resize(0);

                        for (int tr = 0; tr < allTriangles[o].size(); tr++) {

                            glm::vec4 v0 = allTriangles[o][tr].p1();
                            glm::vec4 v1 = allTriangles[o][tr].p2();
                            glm::vec4 v2 = allTriangles[o][tr].p3();

                            if (intersectRayTriangleBarycentric(r.origin(), r.direction(), v0, v1, v2, t_ray)) {
                                if (t_ray < oldt_ray) {
                                    oldt_ray = t_ray;

                                    //Change the colors of the faces to what is in mtl file
                                    if (allobjects[o].mtlYesNo() == 1) {
                                        diffuse_color = allTriDiffColor[o][tr];
                                        specular_color = allTriSpecColor[o][tr];
                                    }

                                    withphong = phongIlluminationonTriangle(r, sourcelight, v0, v1, v2, diffuse_color, specular_color,
                                        r.at(t_ray));

                                    //Check to see if there is any object covering the path from lightsource to camerascene intersection
                                    withphong = shadowTest(r, sourcelight, allTriangles[o], v0, v1, v2, t_ray, withphong);
                                    //change shadow test for array vectors

                                    pixel_color = withphong;
                                }
                            }
                        }

                    }
                }

                //std::cout << lower_left_corner << std::endl;

                int ir = static_cast<int>(255.999 * pixel_color.r);
                int ig = static_cast<int>(255.999 * pixel_color.g);
                int ib = static_cast<int>(255.999 * pixel_color.b);

                outfile << ir << ' ' << ig << ' ' << ib << '\n';

            }
        }

        outfile.close();
        std::cerr << "\nDone.\n";
        printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    }
}