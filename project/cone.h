// cone.h
// test
//
// Created by Nazirul Hasan on 4/10/23.

#ifndef cone_h
#define cone_h

#include <glad/glad.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"

using namespace std;

class Cone {
public:

    // Material properties
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;

    // Texture properties
    float TXmin = 0.0f;
    float TXmax = 1.0f;
    float TYmin = 0.0f;
    float TYmax = 1.0f;
    unsigned int diffuseMap;
    unsigned int specularMap;
    int face;
    float edgeColor[3];

    // Common property
    float shininess;

    // Constructors
    Cone() {
        setUpConeVertexDataAndConfigureVertexAttribute();
    }

    Cone(glm::vec3 amb, glm::vec3 diff, glm::vec3 spec, float shiny) {
        this->ambient = amb;
        this->diffuse = diff;
        this->specular = spec;
        this->shininess = shiny;

        setUpConeVertexDataAndConfigureVertexAttribute();
    }

    Cone(unsigned int dMap, unsigned int sMap, float shiny, float textureXmin, float textureYmin, float textureXmax, float textureYmax, int face = 1) {
        this->diffuseMap = dMap;
        this->specularMap = sMap;
        this->shininess = shiny;
        this->TXmin = textureXmin;
        this->TYmin = textureYmin;
        this->TXmax = textureXmax;
        this->TYmax = textureYmax;
        this->face = face;

        setUpConeVertexDataAndConfigureVertexAttribute();
    }

    // Destructor
    ~Cone() {
        glDeleteVertexArrays(1, &coneVAO);
        glDeleteVertexArrays(1, &lightConeVAO);
        glDeleteVertexArrays(1, &lightTexConeVAO);
        glDeleteBuffers(1, &coneVBO);
        glDeleteBuffers(1, &coneEBO);
    }

    void drawConeWithTexture(Shader& lightingShaderWithTexture, glm::mat4 model = glm::mat4(1.0f)) {
        lightingShaderWithTexture.use();

        // Set material properties
        lightingShaderWithTexture.setInt("material.diffuse", 0);
        lightingShaderWithTexture.setInt("material.specular", 1);
        lightingShaderWithTexture.setFloat("material.shininess", this->shininess);

        // Bind diffuse map
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, this->diffuseMap);

        // Bind specular map
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, this->specularMap);

        // Set model transformation matrix
        lightingShaderWithTexture.setMat4("model", model);

        // Bind the VAO for the textured cone and draw
        glBindVertexArray(lightTexConeVAO);
        glDrawElements(GL_TRIANGLES, 10, GL_UNSIGNED_INT, 0);
    }

    void drawConeWithMaterialisticProperty(Shader& lightingShader, glm::mat4 model = glm::mat4(1.0f)) {
        lightingShader.use();

        // Set material properties
        lightingShader.setVec3("material.ambient", this->ambient);
        lightingShader.setVec3("material.diffuse", this->diffuse);
        lightingShader.setVec3("material.specular", this->specular);
        lightingShader.setFloat("material.shininess", this->shininess);

        // Set model transformation matrix
        lightingShader.setMat4("model", model);

        // Bind the VAO for the cone and draw
        glBindVertexArray(lightConeVAO);
        glDrawElements(GL_TRIANGLES, 10, GL_UNSIGNED_INT, 0);
    }

    void drawCone(Shader& shader, glm::mat4 model = glm::mat4(1.0f), float r = 1.0f, float g = 1.0f, float b = 1.0f) {
        shader.use();

        // Set uniform color for the cone
        shader.setVec3("color", glm::vec3(r, g, b));

        // Set model transformation matrix
        shader.setMat4("model", model);

        // Bind the VAO for the cone and draw
        glBindVertexArray(coneVAO);
        glDrawElements(GL_TRIANGLES, 10, GL_UNSIGNED_INT, 0);
    }

    void setMaterialisticProperty(glm::vec3 amb, glm::vec3 diff, glm::vec3 spec, float shiny) {
        this->ambient = amb;
        this->diffuse = diff;
        this->specular = spec;
        this->shininess = shiny;
    }

    void setTextureProperty(unsigned int dMap, unsigned int sMap, float shiny) {
        this->diffuseMap = dMap;
        this->specularMap = sMap;
        this->shininess = shiny;
    }

private:
    unsigned int coneVAO;
    unsigned int lightConeVAO;
    unsigned int lightTexConeVAO;
    unsigned int coneVBO;
    unsigned int coneEBO;

    void setUpConeVertexDataAndConfigureVertexAttribute() {
        float cone_vertices[] = {
            // Base circle vertices (positions, normals, texture coordinates)
            0.0f, 0.0f, 0.0f,  0.0f, -1.0f, 0.0f,  0.5f, 0.5f, // Center of the base
            1.0f, 0.0f, 0.0f,  0.0f, -1.0f, 0.0f,  1.0f, 0.5f, // First point on the circumference
            0.707f, 0.707f, 0.0f,  0.0f, -1.0f, 0.0f,  0.85f, 0.85f, // Second point on the circumference
            -0.707f, 0.707f, 0.0f,  0.0f, -1.0f, 0.0f,  0.15f, 0.85f, // Third point on the circumference
            -1.0f, 0.0f, 0.0f,  0.0f, -1.0f, 0.0f,  0.0f, 0.5f, // Fourth point on the circumference
            -0.707f, -0.707f, 0.0f,  0.0f, -1.0f, 0.0f,  0.15f, 0.15f, // Fifth point on the circumference
            0.707f, -0.707f, 0.0f,  0.0f, -1.0f, 0.0f,  0.85f, 0.15f, // Sixth point on the circumference
            1.0f, 0.0f, 0.0f,  0.0f, -1.0f, 0.0f,  1.0f, 0.5f, // Back to the first point to close the circle

            // Apex of the cone
            0.0f, 0.0f, 1.0f,  0.0f, 0.0f, 1.0f,  0.5f, 0.5f // Apex
        };

        unsigned int cone_indices[] = {
            // Base (forming a triangle fan with the center point)
            0, 1, 2,
            0, 2, 3,
            0, 3, 4,
            0, 4, 5,
            0, 5, 6,
            0, 6, 1,

            // Side faces (forming triangles between consecutive base vertices and the apex)
            1, 2, 7,  // First side
            2, 3, 7,  // Second side
            3, 4, 7,  // Third side
            4, 5, 7,  // Fourth side
            5, 6, 7,  // Fifth side
            6, 1, 7   // Sixth side
        };


        glGenVertexArrays(1, &coneVAO);
        glGenVertexArrays(1, &lightConeVAO);
        glGenVertexArrays(1, &lightTexConeVAO);
        glGenBuffers(1, &coneVBO);
        glGenBuffers(1, &coneEBO);


        glBindVertexArray(lightTexConeVAO);

        glBindBuffer(GL_ARRAY_BUFFER, coneVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cone_vertices), cone_vertices, GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coneEBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cone_indices), cone_indices, GL_STATIC_DRAW);

        // position attribute
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        // vertex normal attribute
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)12);
        glEnableVertexAttribArray(1);

        // texture coordinate attribute
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)24);
        glEnableVertexAttribArray(2);


        glBindVertexArray(lightConeVAO);

        glBindBuffer(GL_ARRAY_BUFFER, coneVBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coneEBO);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)12);
        glEnableVertexAttribArray(1);


        glBindVertexArray(coneVAO);

        glBindBuffer(GL_ARRAY_BUFFER, coneVBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coneEBO);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }

};


#endif /* cube_h */
#pragma once
