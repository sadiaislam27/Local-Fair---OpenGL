#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define STB_IMAGE_IMPLEMENTATION
#pragma warning(disable:4996)
#include <unordered_map>
#include "shader.h"
#include "camera.h"
#include "basic_camera.h"
#include "pointLight.h"
#include "cube.h"
#include "cone.h"
#include "sphere.h"
#include "stb_image.h"
#include <iostream>
#include <ctime>

using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model, float r, float g, float b);
void getCurrentTime(int& hours, int& minutes, int& seconds);
glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz);

void drawCone(Shader& lightingShader, glm::vec3 color, glm::vec3 translation, glm::vec3 scale, glm::vec3 rot, float radius, float height);
void shop1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& shop1roof, Cube& shop1stick, Cube& shop1base);
void shop2(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& shop1roof, Cube& shop1stick, Cube& shop1base);
void shop(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& shop1roof, Cube& shop1stick, Cube& shop1base, Cube& pastry);
//void dfence(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cone& cone,glm::vec4 color);
void f_dining_room_wall(unsigned int& cubeVAO, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube,glm::vec3 color);
// ************************DRAWING ROOM****************************************
void fan(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& b, Cube& c, float x = 0.0f, float y = 0.0f, float z = 0.0f);
void stage(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Cube& stagewall);
void f_dining_jug(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, float tx, float ty, float tz, float sx, float sy, float sz);

//  *************************************KITCHEN***********************************************
void dfloor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& floor);
void gate(unsigned int& cubeVAO, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Cube& cube1, glm::vec3 color1, glm::vec3 color2);

void drawCylinder(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether);
void canopi(float rotz, Shader& lightingShader, glm::vec3 color, glm::mat4 altogether);
void canopistall(Shader& lightingShaderWithTexture, glm::vec3 color, glm::mat4 alTogether, unsigned int& cubeVAO);
void callcanopistall(Shader& lightingShaderWithTexture, glm::vec3 color, glm::mat4 alTogether, unsigned int& cubeVAO);
void canopi2(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether, GLuint textureID, GLuint textureID2);

void jump(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);

void tree(Shader& lightingShaderWithTexture, glm::vec3 color, glm::mat4 alTogether, unsigned int& cubeVAO, GLuint textureID, GLuint textureID2);


void copyShop1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase);
void copyShop(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2, Cube& juice, Cube& roof1, Cube& roof2);

void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);
void drawLine(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2);
void drawFlywheel(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);
void drawtrain(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);
void drawSeat(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);
void drawCircleflat(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);

void load_texture(unsigned int& texture, string image_name, GLenum format);
unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax);
void shaderActivate(Shader& shader);

glm::mat4 transform(float tr_x, float tr_y, float tr_z, float rot_x, float rot_y, float rot_z, float scal_x, float scal_y, float scal_z) {
    // Modelling Transformation
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tr_x, tr_y, tr_z));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rot_x), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rot_y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rot_z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(scal_x, scal_y, scal_z));
    model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
    return model;
}

// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 800;

// modelling transform
float rotateAngle_X = 0.0;
float rotateAngle_Y = 0.0;
float rotateAngle_Z = 0.0;
float rotateAxis_X = 0.0;
float rotateAxis_Y = 0.0;
float rotateAxis_Z = 1.0;
float translate_X = 0.0;
float translate_Y = 0.0;
float translate_Z = 0.0;
float scale_X = 1.0;
float scale_Y = 1.0;
float scale_Z = 1.0;

// front door
float f_door = 0.0f;
float b1_door = 0.0f;
// bathroom door
float max_bathroom_door_translate = (3.0 + 0.1 - 1.75) / 2.0;
float bathroom_door_translate = 0.0f;

// fan
float rotateFan = 0;
float rotateClock = 0.0f;
bool sign = 1;
bool fanSwitch = 1;

//tvchannel
int channel = 0;

unsigned int sq_tile_tex;
unsigned int black_tex;
unsigned int ch_wood_tex;
unsigned int leaf_tex;


// camera
Camera camera(glm::vec3(0.0f, 1.1f, 5.2f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float eyeX = 0.0, eyeY = 1.0, eyeZ = 3.0;
float lookAtX = 0.0, lookAtY = 0.0, lookAtZ = 0.0;
glm::vec3 V = glm::vec3(0.0f, 1.0f, 0.0f);
BasicCamera basic_camera(eyeX, eyeY, eyeZ, lookAtX, lookAtY, lookAtZ, V);


// positions of the point lights
glm::vec3 pointLightPositions[] = {
    glm::vec3(0.17f, 0.4f, -1.75f),
    glm::vec3(0.0f,  1.5f,  0.0f),
    glm::vec3(0.0f,  1000.0f,  0.0f),
    glm::vec3(0.0f,  3.0f,  0.0f)
};

glm::vec3 point_light_positions[] = {
    glm::vec3(1.45f, 1.3f, 0.1f),
    glm::vec3(1.45f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, 0.5f),
    glm::vec3(2.5f + 1.9f, 0.8f, -0.9f)
};

vector<float>baloon_shape_vertices = {
    -0.0550, 2.3550, 5.1000,
    -1.1900, 1.3050, 5.1000,
    -0.0500, 0.1700, 5.1000,
};

PointLight pointlight1(

    pointLightPositions[0].x, pointLightPositions[0].y, pointLightPositions[0].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    1       // light number
);
PointLight pointlight2(

    pointLightPositions[1].x, pointLightPositions[1].y, pointLightPositions[1].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    2       // light number
);

PointLight pointlight3(

    pointLightPositions[2].x, pointLightPositions[2].y, pointLightPositions[2].z,  // position
    0.1f, 0.1f, 0.1f,     // ambient
    0.1f, 0.1f, 0.1f,      // diffuse
    0.1f, 0.1f, 0.1f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    3       // light number
);
PointLight pointlight4(

    pointLightPositions[3].x, pointLightPositions[3].y, pointLightPositions[3].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    4       // light number
);
// ******************************DRAWING_ROOM_LIGHT***********************************
PointLight drawing_light(
    point_light_positions[0].x, point_light_positions[0].y, point_light_positions[0].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    5       // light number
);
PointLight bed_room1_light(
    point_light_positions[1].x, point_light_positions[1].y, point_light_positions[1].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    6       // light number
);
PointLight dining_light(
    point_light_positions[2].x, point_light_positions[2].y, point_light_positions[2].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    7       // light number
);
PointLight bed_room2_light(
    point_light_positions[3].x, point_light_positions[3].y, point_light_positions[3].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    8       // light number
);
PointLight bathroom_light(
    point_light_positions[4].x, point_light_positions[4].y, point_light_positions[4].z,  // position
    0.2f, 0.2f, 0.2f,     // ambient
    0.2f, 0.2f, 0.2f,      // diffuse
    0.2f, 0.2f, 0.2f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    9       // light number
);



// light settings
bool onOffPointToggle = true;
bool onOffSpotToggle = false;
bool onOffDirectToggle = false;
bool ambientToggle = true;
bool diffuseToggle = true;
bool specularToggle = true;

bool isFlywheelRotating = false;
float flywheelRotationSpeed = 10.0f;
float flywheelRotation = 1;




//glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
//glm::mat4 view = camera.GetViewMatrix();
glm::mat4 projection;
glm::mat4 view;

string diffuseMapPath;
string specularMapPath;
GLuint textureID;
GLuint textureID2;
GLuint textureID3;
GLuint textureID4;



class Curve
{
public:
    vector<float> cntrlPoints;
    vector <float> coordinates;
    vector <float> normals;
    vector <int> indices;
    vector <float> vertices;
    vector<float> texCoords;
    const double pi = 3.14159265389;
    const int nt = 40;
    const int ntheta = 20;
    // Texture properties
    unsigned int diffuseMap;
    unsigned int specularMap;
    float shininess;
    Curve(vector<float>& tmp, unsigned int dMap, unsigned int sMap, float shiny)
        : diffuseMap(dMap), specularMap(sMap), shininess(shiny)
    {
        this->cntrlPoints = tmp;
        this->fishVAO = hollowBezier(cntrlPoints.data(), ((unsigned int)cntrlPoints.size() / 3) - 1);
        cout << cntrlPoints.size() << endl;
        cout << coordinates.size() << endl;
        cout << normals.size() << endl;
        cout << indices.size() << endl;
        cout << vertices.size() << endl;
    }
    ~Curve()
    {
        glDeleteVertexArrays(1, &fishVAO);
        glDeleteVertexArrays(1, &bezierVAO);
        glDeleteBuffers(1, &bezierVBO);
        glDeleteBuffers(1, &bezierEBO);
    }
    void draw(Shader& lightingShader, glm::mat4 model, glm::vec3 amb = glm::vec3(1.0f, 1.0f, 1.0f))
    {
        lightingShader.use();
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("material.ambient", amb);
        lightingShader.setVec3("material.diffuse", amb);
        lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
        lightingShader.setFloat("material.shininess", 32.0f);

        // Set texture properties
        lightingShader.setInt("material.diffuseMap", 0);  // 0 corresponds to GL_TEXTURE0
        lightingShader.setInt("material.specularMap", 1); // 1 corresponds to GL_TEXTURE1

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, diffuseMap);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, specularMap);

        glBindVertexArray(fishVAO);
        glDrawElements(GL_TRIANGLES, (unsigned int)indices.size(), GL_UNSIGNED_INT, (void*)0);

        // unbind VAO
        glBindVertexArray(0);
    }
    void setTextureProperty(unsigned int dMap, unsigned int sMap, float shiny)
    {
        this->diffuseMap = dMap;
        this->specularMap = sMap;
        this->shininess = shiny;
    }


private:
    unsigned int fishVAO;
    unsigned int bezierVAO;
    unsigned int bezierVBO;
    unsigned int bezierEBO;


    unsigned int drawControlPoints()
    {
        unsigned int controlPointVAO;
        unsigned int controlPointVBO;

        glGenVertexArrays(1, &controlPointVAO);
        glGenBuffers(1, &controlPointVBO);

        glBindVertexArray(controlPointVAO);

        glBindBuffer(GL_ARRAY_BUFFER, controlPointVBO);
        glBufferData(GL_ARRAY_BUFFER, (unsigned int)cntrlPoints.size() * sizeof(float), cntrlPoints.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        return controlPointVAO;
    }

    long long nCr(int n, int r)
    {
        if (r > n / 2)
            r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for (i = 1; i <= r; i++)
        {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }
    void BezierCurve(double t, float xy[2], GLfloat ctrlpoints[], int L)
    {
        double y = 0;
        double x = 0;
        t = t > 1.0 ? 1.0 : t;
        for (int i = 0; i < L + 1; i++)
        {
            long long ncr = nCr(L, i);
            double oneMinusTpow = pow(1 - t, double(L - i));
            double tPow = pow(t, double(i));
            double coef = oneMinusTpow * tPow * ncr;
            x += coef * ctrlpoints[i * 3];
            y += coef * ctrlpoints[(i * 3) + 1];

        }
        xy[0] = float(x);
        xy[1] = float(y);
    }
    unsigned int hollowBezier(GLfloat ctrlpoints[], int L)
    {
        int i, j;
        float x, y, z, r;                // current coordinates
        float theta;
        float nx, ny, nz, lengthInv;    // vertex normal
        float u, v;                     // texture coordinates

        const float dtheta = 2 * pi / ntheta;        // angular step size

        float t = 0;
        float dt = 1.0 / nt;
        float xy[2];

        for (i = 0; i <= nt; ++i)              // step through y
        {
            BezierCurve(t, xy, ctrlpoints, L);
            r = xy[0];
            y = xy[1];
            theta = 0;
            t += dt;
            lengthInv = 1.0 / r;

            for (j = 0; j <= ntheta; ++j)
            {
                double cosa = cos(theta);
                double sina = sin(theta);
                z = r * cosa;
                x = r * sina;

                coordinates.push_back(x);
                coordinates.push_back(y);
                coordinates.push_back(z);

                // normalized vertex normal (nx, ny, nz)
                // center point of the circle (0,y,0)
                nx = (x - 0) * lengthInv;
                ny = (y - y) * lengthInv;
                nz = (z - 0) * lengthInv;

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);

                // Calculate texture coordinates (s, t)
                u = static_cast<float>(j) / ntheta;
                v = static_cast<float>(i) / nt;
                texCoords.push_back(u);
                texCoords.push_back(v);

                theta += dtheta;
            }
        }
        // generate index list of triangles
        // k1--k1+1
        // |  / |
        // | /  |
        // k2--k2+1

        int k1, k2;
        for (int i = 0; i < nt; ++i)
        {
            k1 = i * (ntheta + 1);     // beginning of current stack
            k2 = k1 + ntheta + 1;      // beginning of next stack

            for (int j = 0; j < ntheta; ++j, ++k1, ++k2)
            {
                // k1 => k2 => k1+1
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);

                // k1+1 => k2 => k2+1
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }

        size_t count = coordinates.size();
        for (i = 0, j = 0; i < count; i += 3, j += 2)
        {
            //cout << count << ' ' << i + 2 << endl;
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            if (i < normals.size())
                vertices.push_back(normals[i]);
            if (i + 1 < normals.size())
                vertices.push_back(normals[i + 1]);
            if (i + 2 < normals.size())
                vertices.push_back(normals[i + 2]);

            //// Add texture coordinates
            //if (j < texCoords.size())
            //    vertices.push_back(texCoords[j]);
            //if (j + 1 < texCoords.size())
            //    vertices.push_back(texCoords[j + 1]);
        }

        glGenVertexArrays(1, &bezierVAO);
        glBindVertexArray(bezierVAO);

        // create VBO to copy vertex data to VBO
        glGenBuffers(1, &bezierVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bezierVBO);           // for vertex data
        glBufferData(GL_ARRAY_BUFFER,                   // target
            (unsigned int)vertices.size() * sizeof(float), // data size, # of bytes
            vertices.data(),   // ptr to vertex data
            GL_STATIC_DRAW);                   // usage

        // create EBO to copy index data
        glGenBuffers(1, &bezierEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bezierEBO);   // for index data
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,           // target
            (unsigned int)indices.size() * sizeof(unsigned int),             // data size, # of bytes
            indices.data(),               // ptr to index data
            GL_STATIC_DRAW);                   // usage

        // activate attrib arrays
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
        glEnableVertexAttribArray(2);
        // set attrib arrays with stride and offset
        int stride = 24;  // should be 24 bytes
        glVertexAttribPointer(0, 3, GL_FLOAT, false, stride, (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, false, stride, (void*)(sizeof(float) * 3));
        glVertexAttribPointer(2, 2, GL_FLOAT, false, stride, (void*)(sizeof(float) * 6)); // Add this line for texture coordinates
        //
                // unbind VAO, VBO and EBO
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return bezierVAO;
    }

};
Curve* bucket;
Curve* leaf;
Curve* wood;

vector<float>Bucket = {
0.001, 1.7054, 4.9697,
-0.550, 0.5803, 5.2077,
-1.00, -0.0164, 5.3340,
-1.500, -0.2268, 5.3785,

};

vector<float>tree2_vertices = {
    0.0250, 1.6600, 5.1000,
    -0.250, 1.3900, 5.1000,
    -0.3450, 1.3550, 5.1000,
    -0.3850, 1.2100, 5.1000,
    -1.30, 0.6700, 5.1000,
    -.850, .100, 5.1000,
    -1.20, 0.200, 5.1000,
    -0.3250, 0.2300, 5.1000
};
/*vector<float>tree2_vertices = {
    0.0250, 1.6600, 5.1000,
    -0.7250, 1.3900, 5.1000,
    -0.2450, 1.3550, 5.1000,
    -0.8500, 1.2300, 5.1000,
    -0.2850, 1.2100, 5.1000,
    -0.8650, 0.7700, 5.1000,
    -0.3250, 0.8300, 5.1000
};*/



// timing
float deltaTime = 0.0f;    // time between current frame and last frame
float lastFrame = 0.0f;

unsigned int texture0, texture1;

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Local Fair", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader zprogram
    // ------------------------------------
    Shader lightingShader("vertexShaderForPhongShading.vs", "fragmentShaderForPhongShading.fs");
    Shader lightingShaderWithTexture("vertexShaderForPhongShadingWithTexture.vs", "fragmentShaderForPhongShadingWithTexture.fs");
    Shader ourShader("vertexShader.vs", "fragmentShader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------

    float cube_vertices[] = {
        // positions      // normals
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,

        1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,

        0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,

        1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,

        0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f
    };
    unsigned int cube_indices[] = {
        0, 3, 2,
        2, 1, 0,

        4, 5, 7,
        7, 6, 4,

        8, 9, 10,
        10, 11, 8,

        12, 13, 14,
        14, 15, 12,

        16, 17, 18,
        18, 19, 16,

        20, 21, 22,
        22, 23, 20
    };

    unsigned int cubeVAO, cubeVBO, cubeEBO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vertices), cube_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
    unsigned int lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glBindVertexArray(lightCubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    // note that we update the lamp's position attribute's stride to reflect the updated buffer data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    //stall1roof

    


    diffuseMapPath = "images/roof_shop.jpg";
    specularMapPath = "images/roof_shop.jpg";
    Cube stall1roof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);
    
    diffuseMapPath = "images/roof1.jpg";
    specularMapPath = "images/roof1.jpg";
    Cube roof1 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/roof9.jpg";
    specularMapPath = "images/roof9.jpg";
    Cube roof2 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/roof8.jpg";
    specularMapPath = "images/roof8.jpg";
    Cube roof3 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/roof7.jpg";
    specularMapPath = "images/roof7.jpg";
    Cube roof4 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/roof4.jpg";
    specularMapPath = "images/roof4.jpg";
    Cube roof5 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/canopi1.png";
    black_tex = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    Curve buck(Bucket, black_tex, black_tex, 1.0f);
    bucket = &buck;


    diffuseMapPath = "images/tree.jpg";
    leaf_tex = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Curve Leaf(tree2_vertices, leaf_tex, leaf_tex, 1.0f);
    leaf = &Leaf;

    //baloon_shape_vertices

    textureID3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    textureID4 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    diffuseMapPath = "images/wood2.jpg";
    ch_wood_tex = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Curve Wood(baloon_shape_vertices, ch_wood_tex, ch_wood_tex, 1.0f);
    wood = &Wood;

    texture0 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    texture1 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    diffuseMapPath = "images/roof_shop1.png";
    specularMapPath = "images/roof_shop1.png";
    Cube stall2roof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    textureID = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    textureID2 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    diffuseMapPath = "images/stick1.jpg";
    specularMapPath = "images/stick1.jpg";
    Cube stallstick = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);


    diffuseMapPath = "images/wood2.jpg";
    specularMapPath = "images/wood2.jpg";
    Cube stallbase = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/grass2.jpg";
    specularMapPath = "images/grass2.jpg";
    Cube floor = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);


    diffuseMapPath = "images/fence4.jpg";
    specularMapPath = "images/fence4.jpg";
    Cube fence = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 10.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/curtain.png";
    specularMapPath = "images/curtain.png";
    Cube stagewall = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/wall_color.png";
    specularMapPath = "images/wall_color.png";
    Cube stagew = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);
    
    diffuseMapPath = "images/wood2.jpg";
    specularMapPath = "images/wood2.jpg";
    Cone cone = Cone(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, .0f, 6.0f, 4.0f);


    diffuseMapPath = "images/gate.jpg";
    specularMapPath = "images/gate.jpg";
    Cube Gate = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);
    

    diffuseMapPath = "images/welcome3.png";
    specularMapPath = "images/welcome3.png";
    Cube welcome  = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/pastry2.jpg";
    specularMapPath = "images/pastry2.jpg";
    Cube pastry = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/juice.jpg";
    specularMapPath = "images/juice.jpg";
    Cube juice = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    
    
    //ourShader.use();
    //lightingShader.use();

    

    pointlight2.turnOff();

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.5f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShader.use();
        lightingShader.setVec3("viewPos", camera.Position);

        // point light 1
        pointlight1.setUpPointLight(lightingShader);
        // point light 2
        pointlight2.setUpPointLight(lightingShader);

        lightingShader.setVec3("diectionalLight.directiaon", 0.0f, 3.0f, 0.0f);
        lightingShader.setVec3("diectionalLight.ambient", .2, .2, .2);
        lightingShader.setVec3("diectionalLight.diffuse", .8f, .8f, .8f);
        lightingShader.setVec3("diectionalLight.specular", 1.0f, 1.0f, 1.0f);

        lightingShader.setBool("dlighton", false);


        lightingShader.setVec3("spotlight.position", 0.5, 1, -0.5);
        lightingShader.setVec3("spotlight.direction", 0, -1, 0);
        lightingShader.setVec3("spotlight.ambient", .2, .2, .2);
        lightingShader.setVec3("spotlight.diffuse", .8f, .8f, .8f);
        lightingShader.setVec3("spotlight.specular", 1.0f, 1.0f, 1.0f);
        lightingShader.setFloat("spotlight.k_c", 1.0f);
        lightingShader.setFloat("spotlight.k_l", 0.09);
        lightingShader.setFloat("spotlight.k_q", 0.032);
        lightingShader.setFloat("cos_theta", glm::cos(glm::radians(5.5f)));
        lightingShader.setBool("spotlighton", false);

        // activate shader
        lightingShader.use();

        projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        lightingShader.setMat4("projection", projection);

        // camera/view transformation
        // glm::mat4 view = camera.GetViewMatrix();
        //glm::mat4 view = basic_camera.createViewMatrix();
        view = camera.GetViewMatrix();
        lightingShader.setMat4("view", view);

        // Modelling Transformation
        glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShader.setMat4("model", model);

        float r = 0.5, h = 0.5;
        //drawCone(lightingShader, glm::vec3(0.1, 0.2, 0.2), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(rotateAngle_X, rotateAngle_Y, rotateAngle_Z), r, h);

        if (!fanSwitch)rotateFan += 35.0f;

        //DRAWING ROOM

        //  drawing_tv
        //cout << channel << " ";
        if(channel == 0) {
            diffuseMapPath = "images/channel3.jpg";
            specularMapPath = "images/channel3.jpg";
        }
        else if (channel == 1) {
            diffuseMapPath = "images/channel1.jpg";
            specularMapPath = "images/channel1.jpg";
        }
        else if (channel == 2) {
            diffuseMapPath = "images/channel6.jpg";
            specularMapPath = "images/channel6.jpg";
        }
        else if (channel == 3) {
            diffuseMapPath = "images/channel5.jpg";
            specularMapPath = "images/channel5.jpg";
        }
        
        Cube drawing_tv = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


        glm::vec3 color = glm::vec3(1.0, 1.0,1.0);
        f_dining_room_wall(cubeVAO,  model, lightingShaderWithTexture, fence,color);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

        glm::vec3 color1 = glm::vec3(1.0, 1.0, 1.0);
        glm::vec3 color2 = glm::vec3(1.0, 1.0, 1.0);
        //gate(cubeVAO, model, lightingShaderWithTexture, Gate, welcome, color1,color2);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        //canopistall(lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f), model,cubeVAO);
        callcanopistall(lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f), model, cubeVAO);

        //tree(lightingShaderWithTexture, glm::vec3(.0f, .5f, .0f), model,cubeVAO, texture0, texture1);

        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

       // canopi2(lightingShaderWithTexture, glm::vec3(1.0f, 1.0f, 1.0f), model, textureID, textureID2);
        
        
        //stage(cubeVAO, lightingShader, model, lightingShaderWithTexture, stagew,stagewall);
        dfloor(cubeVAO, lightingShader, model, lightingShaderWithTexture, floor);

        //drawFlywheel(model, model, lightingShaderWithTexture, cubeVAO);
        
        //shop(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick,stallbase, pastry);
        //shop1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall2roof, stagew, stallbase);
        
        copyShop(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, pastry,juice,roof1,roof2);
        copyShop1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall2roof, stagew, stallbase);

        glm::vec3 flywheelPosition(translate_X, translate_Y, translate_Z);
        glm::vec3 flywheelScale(scale_X, scale_Y, scale_Z);     // Scale

        glm::vec3 flywheelColor(.3f, .50f, .60f);

        glm::mat4 model1, translateM;
        translateM = glm::mat4(1.0);
        translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        model = rotateXMatrix * rotateYMatrix * rotateZMatrix;
        model1 = translateM * rotateXMatrix * rotateYMatrix * rotateZMatrix;
        lightingShaderWithTexture.setMat4("model", model1);
        lightingShaderWithTexture.setVec3("curveColor", flywheelColor);

        drawtrain(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);
        jump(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);
        //drawSeat(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);
        //void jump(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO)
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

        if (isFlywheelRotating)
        {
            glm::vec3 flywheelPosition(translate_X, translate_Y, translate_Z);
            glm::vec3 flywheelScale(scale_X, scale_Y, scale_Z);

            glm::vec3 flywheelColor(0.0f, 0.0f, 10.0f);
            glm::mat4 model1, translateM;
            translateM = glm::mat4(1.0);
            translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
            rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
            rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
            rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
            // glm::mat4 rotateFlywheelMatrix = glm::rotate(identityMatrix, glm::radians(flywheelRotation), glm::vec3(0.0f, 0.0f, 1.0f));
            model = rotateXMatrix * rotateYMatrix * rotateZMatrix;
            model1 = translateM * rotateXMatrix * rotateYMatrix * rotateZMatrix;
            lightingShaderWithTexture.setMat4("model", model);


            drawFlywheel(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);


            flywheelRotation += 1;
        }

        else
        {
            //glm::vec3 flywheelPosition(-4.0f, 2.0f, -4.0f); // Position
            glm::vec3 flywheelPosition(translate_X, translate_Y, translate_Z);
            glm::vec3 flywheelScale(scale_X, scale_Y, scale_Z);     // Scale

            glm::vec3 flywheelColor(0.0f, 0.0f, 10.0f);

            glm::mat4 model1, translateM;
            translateM = glm::mat4(1.0);
            translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
            rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
            rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
            rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
            model1 = translateM * rotateXMatrix * rotateYMatrix * rotateZMatrix;
            lightingShaderWithTexture.setMat4("model", model);


            drawFlywheel(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);



        }


        //model = transform(-24.0f, 1, -26, 0.0f, 0.0f, 0.0f, 1.2, 1.2, 1.2);
        //make_tree2(cube, wheel, wheel_hollow, tree2, lightingShader, lightingShaderWithTexture, model);
        // 
        //the lamp objects
        ourShader.use();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);
        glBindVertexArray(lightCubeVAO);
        

        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X + 1, translate_Y + 1, translate_Z + 1));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X * 0.5, scale_Y * 0.5, scale_Z * 0.5));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

        lightingShaderWithTexture.use();
        lightingShaderWithTexture.setVec3("viewPos", camera.Position);
        lightingShaderWithTexture.setMat4("view", view);
        lightingShaderWithTexture.setMat4("projection", projection);

        lightingShaderWithTexture.use();
        // point light 1
        pointlight1.setUpPointLight(lightingShaderWithTexture);
        // point light 2
        pointlight2.setUpPointLight(lightingShaderWithTexture);
        // point light 3
        pointlight3.setUpPointLight(lightingShaderWithTexture);
        // point light 4
        pointlight4.setUpPointLight(lightingShaderWithTexture);

        drawing_light.setUpPointLight(lightingShaderWithTexture);
        bed_room1_light.setUpPointLight(lightingShaderWithTexture);
        dining_light.setUpPointLight(lightingShaderWithTexture);
        bed_room2_light.setUpPointLight(lightingShaderWithTexture);
        bathroom_light.setUpPointLight(lightingShaderWithTexture);

      
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &cubeVAO);
    glDeleteVertexArrays(1, &lightCubeVAO);
    glDeleteBuffers(1, &cubeVBO);
    glDeleteBuffers(1, &cubeEBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

void shaderActivate(Shader& shader)
{
    shader.use();
    shader.setVec3("viewPos", camera.Position);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);
}


void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model = glm::mat4(1.0f), float r = 1.0f, float g = 1.0f, float b = 1.0f)
{
    lightingShader.use();

    lightingShader.setVec3("material.ambient", glm::vec3(r, g, b));
    lightingShader.setVec3("material.diffuse", glm::vec3(r, g, b));
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);

    lightingShader.setMat4("model", model);

    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides) {
    // Generate vertex data for the circle
    GLint numberOfVertices = numberOfSides + 1;
    GLfloat doublePi = 2.0f * 3.14159265359f;

    std::vector<GLfloat> vertices;

    for (int i = 0; i < numberOfVertices; i++) {
        vertices.push_back(x + (radius * cos(i * doublePi / numberOfSides))); // X coordinate
        vertices.push_back(z + (radius * sin(i * doublePi / numberOfSides))); // Y coordinate
        vertices.push_back(y);                                               // Z coordinate
    }

    // Create and bind a VAO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Bind and upload vertex data to VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    // Define vertex attribute layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Draw the circle as a line loop
    glDrawArrays(GL_LINE_LOOP, 0, numberOfVertices);

    // Cleanup
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}


void drawCircleflat(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides) {
    // Generate vertex data for the circle
    GLint numberOfVertices = numberOfSides + 1;
    GLfloat doublePi = 2.0f * 3.14159265359f;

    std::vector<GLfloat> vertices;

    for (int i = 0; i < numberOfVertices; i++) {

        vertices.push_back(x + (radius * cos(i * doublePi / numberOfSides))); // X coordinate
        vertices.push_back(y); // Y coordinate
        vertices.push_back(z + (radius * sin(i * doublePi / numberOfSides)));                                               // Z coordinate
    }

    // Create and bind a VAO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Bind and upload vertex data to VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    // Define vertex attribute layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Draw the circle as a line loop
    glDrawArrays(GL_LINE_LOOP, 0, numberOfVertices);

    // Cleanup
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}


void drawLine(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2) {
    // Vertex data for the line
    GLfloat vertices[] = {
        x1, y1, z1,
        x2, y2, z2
    };

    // Create and bind a VAO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Bind and upload vertex data to VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // Define vertex attribute layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Draw the line
    glDrawArrays(GL_LINES, 0, 2);

    // Cleanup
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}

void drawCylinder(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether)
{
    lightingShader.use();
    // building model matrix
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, model;
    model = glm::mat4(1.0);
    model = model * altogether;
    lightingShader.setMat4("model", model);
    lightingShader.setVec3("curveColor", color);

    int n = 36;

    vector<float>v;
    vector<int>in;

    for (int i = 0; i < n; i++)
    {
        float ang = (i * 2 * 3.1416) / n;
        float x1 = .5 * cos(ang), z1 = .5 * sin(ang), y1 = .5;
        v.push_back(x1); v.push_back(y1); v.push_back(z1);
        v.push_back(color.r); v.push_back(color.g); v.push_back(color.b);
        float x2 = .5 * cos(ang), z2 = .5 * sin(ang), y2 = -.5;

        v.push_back(x2); v.push_back(y2); v.push_back(z2);
        v.push_back(color.r); v.push_back(color.g); v.push_back(color.b);
    }


    for (int i = 0; i < n - 1; i++)
    {
        in.push_back(i * 2);
        in.push_back((i * 2) + 1);
        in.push_back((i + 1) * 2);
        in.push_back((i + 1) * 2);
        in.push_back((i * 2) + 1);
        in.push_back((i + 1) * 2 + 1);

    }

    in.push_back((n - 1) * 2);
    in.push_back((n - 1) * 2 + 1);
    in.push_back(0);
    in.push_back(0);
    in.push_back((n - 1) * 2 + 1);
    in.push_back(1);



    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, in.size() * sizeof(float), in.data(), GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // setting up materialistic property
    lightingShader.setVec3("material.ambient", color);
    lightingShader.setVec3("material.diffuse", color);
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);




    glDrawElements(GL_TRIANGLES, in.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteVertexArrays(1, &VAO);
}

void canopi(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether) {

    shaderActivate(lightingShader);
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    model = identityMatrix * altogether;
    lightingShader.setMat4("model", model);

    lightingShader.setVec3("curveColor", color);
    bucket->draw(lightingShader, model);
    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

}

void copyShop1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase)
{
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model, rotateMatrix;
    glm::mat4 identityMatrix = glm::mat4(1.0f);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-1.60, -.0, -6.0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, 1.0, 1.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    shop1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(.20, -.0, -6.0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, 1.0, 1.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    shop1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-4.10, -.0, -6.0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, 1.0, 1.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    //shop1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase);



}

void copyShop(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2, Cube& juice, Cube& roof1, Cube& roof2)
{
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model, rotateMatrix;
    glm::mat4 identityMatrix = glm::mat4(1.0f);



    translateMatrix = glm::translate(identityMatrix, glm::vec3(-.50, -.0, .75));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, 1.0, 1.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    shop(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-.50, -.0, -1.650));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, 1.0,1.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    shop(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, juice);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-.50, -.0, -4.10));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0,1.0,1.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    shop(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);

}

void canopi2(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether, GLuint textureID, GLuint textureID2) {
    // Activate the shader
    shaderActivate(lightingShader);

    // Identity matrix to reset transformations
    glm::mat4 identityMatrix = glm::mat4(1.0f);

    // Apply transformations
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    model = identityMatrix * altogether;
    lightingShader.setMat4("model", model);

    // Set the curve color
    lightingShader.setVec3("curveColor", color);

    // Bind the texture
    glActiveTexture(GL_TEXTURE0);  // Use texture unit 0
    glBindTexture(GL_TEXTURE_2D, textureID);  // Bind the texture

    // Set the texture uniform in the shader
    lightingShader.setInt("material.diffuse", 0); // Use texture unit 0 for diffuse map

    glActiveTexture(GL_TEXTURE0);  // Use texture unit 0
    glBindTexture(GL_TEXTURE_2D, textureID2);  // Bind the texture

    // Set the texture uniform in the shader
    lightingShader.setInt("material.specular", 1);
    
    // Draw the object (bucket)
    leaf->draw(lightingShader, model);

    // Reset curve color to white
    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));
}


void tree(Shader& lightingShaderWithTexture, glm::vec3 color, glm::mat4 alTogether, unsigned int& cubeVAO, GLuint textureID, GLuint textureID2)
{

    shaderActivate(lightingShaderWithTexture);

    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-.5, 0.8, 0.0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0, .4, .6));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);

    // Draw canopi
    lightingShaderWithTexture.setVec3("curveColor", color);
    canopi2(lightingShaderWithTexture, color, model,texture0,texture1);

}

void canopistall(Shader& lightingShaderWithTexture, glm::vec3 color, glm::mat4 alTogether, unsigned int& cubeVAO)
{

    shaderActivate(lightingShaderWithTexture);

    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-3.5, 0.8, -5.0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.6, .4, .6));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);

    // Draw canopi
    lightingShaderWithTexture.setVec3("curveColor", color);
    canopi(lightingShaderWithTexture, glm::vec3(1.0, 0.85, 0.7), model);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-3.56, -0.8, -5.05));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.1, 2.2, .1));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model );
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.8, 1.0, 0.8));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));
    identityMatrix = glm::mat4(1.0f);
    translateMatrix = glm::mat4(1.0f);
    scaleMatrix = glm::mat4(1.0f);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-3.96, -.90, -5.55));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.1, .1, 1.1));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.45, 0.35, 0.25));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

    identityMatrix = glm::mat4(1.0f);
    translateMatrix = glm::mat4(1.0f);
    scaleMatrix = glm::mat4(1.0f);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-3.15, -.90, -5.54));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.3, .6, 1.1));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.902, 0.82, 0.2));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));
}
void callcanopistall(Shader& lightingShaderWithTexture, glm::vec3 color, glm::mat4 alTogether, unsigned int& cubeVAO)
{
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model, rotateMatrix;
    glm::mat4 identityMatrix = glm::mat4(1.0f);



    translateMatrix = glm::translate(identityMatrix, glm::vec3(.20, .30, -.30));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.3, 1.3, 1.3));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    canopistall(lightingShaderWithTexture,color, model, cubeVAO);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-9.10, .450, -3.0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.4, 1.4, 1.4));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-86.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    canopistall(lightingShaderWithTexture, color, model, cubeVAO);
}

void drawFlywheel(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO) {
    // Apply transformations
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 model3 = glm::mat4(1.0f);
    glm::mat4 model4 = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(2.8, 2.1, -1.2));
    model1 = glm::translate(model1, glm::vec3(2.8, 2.1, -1.2));// Apply translation
    // Apply rotation
    model = glm::scale(model, glm::vec3(2, 2, 1));
    model1 = glm::scale(model1, glm::vec3(2, 2, 1));

    model = alTogether * model;
    model1 = alTogether2 * model1;

    if (isFlywheelRotating)
    {
        glm::mat4 rotateFlywheelMatrix = glm::rotate(model4, glm::radians(flywheelRotation), glm::vec3(0.0f, 0.0f, 1.0f));
        model3 = model * rotateFlywheelMatrix;
    }
    else
    {
        model3 = model;
    }
    lightingShaderWithTexture.setMat4("model", model3);

    lightingShaderWithTexture.setVec3("curveColor", color);



    // Set line width

    glLineWidth(10.0f);

    //// Draw circles
    //drawCircle(0.0f, -0.4f, 0.0f, 0.5f, 36);  // Small circle behind
    //drawCircle(0.0f, 0.0f, 0.0f, 0.5f, 36);   // Large circle in front
    //drawCircle(0.0f, -0.4f, 0.0f, 0.9f, 36);  // Larger circle behind
    //drawCircle(0.0f, 0.0f, 0.0f, 0.9f, 36);   // Larger circle in front
    //drawCircle(0.0f, -0.4f, 0.0f, 0.2f, 36);  // Smallest circle behind
    //drawCircle(0.0f, 0.0f, 0.0f, 0.2f, 36);   // Smallest circle in front

    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.6, 0.0, 0.6));
    drawCircle(0.0f, -0.4f, 0.0f, 0.35f, 36);  // Small circle behind
    drawCircle(0.0f, 0.0f, 0.0f, 0.35f, 36);   // Large circle in front
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(.60, 0.0, 0.60));
    drawCircle(0.0f, -0.4f, 0.0f, 0.9f, 36);  // Larger circle behind
    drawCircle(0.0f, 0.0f, 0.0f, 0.9f, 36);   // Larger circle in front
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 0.0));
    drawCircle(0.0f, -0.4f, 0.0f, 0.2f, 36);  // Smallest circle behind
    drawCircle(0.0f, 0.0f, 0.0f, 0.2f, 36);   // Smallest circle in front


    for (int i = 0; i < 20; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 20;
        GLfloat x1 = 0.9f * cos(angleRad), y1 = 0.9f * sin(angleRad), z1 = -0.0f;
        GLfloat x2 = 0.9f * cos(angleRad), y2 = 0.9f * sin(angleRad), z2 = -0.4f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    for (int i = 0; i < 20; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 20;
        GLfloat x1 = 0.35f * cos(angleRad), y1 = 0.35f * sin(angleRad), z1 = -0.0f;
        GLfloat x2 = 0.35f * cos(angleRad), y2 = 0.35f * sin(angleRad), z2 = -0.4f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    // Draw lines connecting smaller circles
    for (int i = 0; i < 6; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 6;
        GLfloat x1 = 0.2f * cos(angleRad), y1 = 0.2f * sin(angleRad), z1 = -0.4f;
        GLfloat x2 = 0.9f * cos(angleRad), y2 = 0.9f * sin(angleRad), z2 = -0.4f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }
    for (int i = 0; i < 6; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 6;
        GLfloat x1 = 0.2f * cos(angleRad), y1 = 0.2f * sin(angleRad), z1 = 0.0f;
        GLfloat x2 = 0.9f * cos(angleRad), y2 = 0.9f * sin(angleRad), z2 = 0.0f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    // Draw lines connecting larger circles
    for (int i = 0; i < 36; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 36;
        GLfloat x1 = 0.0f ;         // Center of the first circle
        GLfloat y1 = 0.0f;
        GLfloat z1 = -0.4f;        // z-position of the first circle
        GLfloat x2 = 0.2f * cos(angleRad);
        GLfloat y2 = 0.2f * sin(angleRad);
        GLfloat z2 = -0.4f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }
    for (int i = 0; i < 36; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 36;
        GLfloat x1 = 0.0f ;         // Center of the second circle
        GLfloat y1 = 0.0f;
        GLfloat z1 = 0.0f;         // z-position of the second circle
        GLfloat x2 = 0.2f * cos(angleRad) ;
        GLfloat y2 = 0.2f * sin(angleRad);
        GLfloat z2 = 0.0f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }
    for (int i = 0; i < 36; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 36;
        GLfloat x1 = 0.2f * cos(angleRad);
        GLfloat y1 = 0.2f * sin(angleRad);
        GLfloat z1 = -0.4f;        // Points on the first circle
        GLfloat x2 = 0.2f * cos(angleRad);
        GLfloat y2 = 0.2f * sin(angleRad);
        GLfloat z2 = 0.0f;         // Corresponding points on the second circle
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    for (int i = 0; i < 6; i++)
    {
        float angle = i * 60.0f; // 60 degrees between each seat
        glm::mat4 seatModel = glm::mat4(1.0f);
        seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
        seatModel = glm::rotate(seatModel, glm::radians(-rotationAngle), glm::vec3(0.0f, 0.0f, 1.0f)); // rotate to match the Ferris wheel's rotation
        seatModel = glm::translate(seatModel, glm::vec3(-0.2f, -0.2f, -0.5f)); // move to the correct height
        seatModel = glm::scale(seatModel, glm::vec3(0.5f, 0.2f, 0.6f));
        lightingShaderWithTexture.setMat4("model", model3 * seatModel);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 0.5, 0.3));
        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        seatModel = glm::mat4(1.0f);
        seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
        seatModel = glm::rotate(seatModel, glm::radians(-rotationAngle), glm::vec3(0.0f, 0.0f, 1.0f)); // rotate to match the Ferris wheel's rotation
        seatModel = glm::translate(seatModel, glm::vec3(-0.2f, -0.2f, -0.5f)); // move to the correct height
        seatModel = glm::scale(seatModel, glm::vec3(0.5f, 0.2f, 0.6f));
        lightingShaderWithTexture.setMat4("model", model3 * seatModel);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 0.8, 0.4));
        glBindVertexArray(cubeVAO);
        //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    }

    glm::mat4 seatModel = glm::mat4(1.0f);
    //seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
    seatModel = glm::rotate(seatModel, glm::radians(20.0f), glm::vec3(1.0f, 0.0f, 0.0f)); // rotate to match the Ferris wheel's rotation
    seatModel = glm::translate(seatModel, glm::vec3(-0.06f, -1.7f, -0.41f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(0.1f, 1.5f, 0.1f));

    lightingShaderWithTexture.setMat4("model", model * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.8, 1.0, 0.8));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


    seatModel = glm::mat4(1.0f);
    //seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
    seatModel = glm::rotate(seatModel, glm::radians(-20.0f), glm::vec3(1.0f, 0.0f, 0.0f)); // rotate to match the Ferris wheel's rotation
    seatModel = glm::translate(seatModel, glm::vec3(-0.06f, -1.5f, -0.05f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(0.1f, 1.5f, 0.1f));

    lightingShaderWithTexture.setMat4("model", model * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.8, 1.0, 0.8));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


    seatModel = glm::mat4(1.0f);
    seatModel = glm::translate(seatModel, glm::vec3(-.30f, -1.51f, -1.35f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(0.55f, .1f, 2.4f));
    lightingShaderWithTexture.setMat4("model", model1 * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.9, 0.6, 0.4));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


}


float cubeYPosition = -0.51f;
float cubeVelocity = 0.0f;
bool isCubeMovingUp = false;
bool isCubeMovingDown = false;
void jump(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO) {
    // Apply transformations
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 model3 = glm::mat4(1.0f);
    glm::mat4 model4 = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(-.20, .1, -6.0));
    model1 = glm::translate(model1, glm::vec3(-.20, .1, -6.0));// Apply translation
    // Apply rotation
    model = glm::scale(model, glm::vec3(2, 2, 1));
    model1 = glm::scale(model1, glm::vec3(2, 2, 1));

    model = alTogether * model;
    model1 = alTogether2 * model1;


    model3 = model;

    lightingShaderWithTexture.setMat4("model", model3);

    lightingShaderWithTexture.setVec3("curveColor", color);

    glLineWidth(10.0f);
    glm::mat4 seatModel = glm::mat4(1.0f);
    seatModel = glm::translate(seatModel, glm::vec3(-.30f, -.51f, -2.1f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(0.4f, .4f, 0.4f));
    lightingShaderWithTexture.setMat4("model", model * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.74, 0.89, 1.0));
    glBindVertexArray(cubeVAO);
    //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0,1.0));
    seatModel = glm::mat4(1.0f);
    seatModel = glm::translate(seatModel, glm::vec3(-.40f, -.51f, -2.35f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(.50f, 3.1f, .2f));
    lightingShaderWithTexture.setMat4("model", model1 * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 0.8));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
    /*
        seatModel = glm::mat4(1.0f);
        seatModel = glm::translate(seatModel, glm::vec3(-.370f, -.51f, -1.6f)); // move to the correct height
        seatModel = glm::scale(seatModel, glm::vec3(.7f, .1f, .7f));
        lightingShaderWithTexture.setMat4("model", model1 * seatModel);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.8, 1.0, 0.8));
        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    */

    if (isCubeMovingUp) {
        cubeVelocity = 0.1f;
        cubeYPosition += cubeVelocity;
        if (cubeYPosition >= 2.0f) { // Change to desired height
            isCubeMovingDown = true;
            isCubeMovingUp = false;
        }
    }
    else if (isCubeMovingDown) {
        cubeVelocity = -0.1f;
        cubeYPosition += cubeVelocity;
        if (cubeYPosition <= -0.1f) {
            isCubeMovingUp = true;
            isCubeMovingDown = false;
        }
    }

    // Draw the cube with the updated position
 
    // Draw the cube with the updated position
    seatModel = glm::mat4(1.0f);
    seatModel = glm::translate(seatModel, glm::vec3(-.490f, cubeYPosition, -2.18f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(.7f, .1f, .7f));
    lightingShaderWithTexture.setMat4("model", model1 * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.1, 0.1, 0.44));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    seatModel = glm::mat4(1.0f);
    seatModel = glm::translate(seatModel, glm::vec3(-.490f, cubeYPosition+0.1, -2.18f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(.1f, .1f, .7f));
    lightingShaderWithTexture.setMat4("model", model1 * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.1, 0.1, 0.44));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    seatModel = glm::mat4(1.0f);
    seatModel = glm::translate(seatModel, glm::vec3(.11f, cubeYPosition + 0.1, -2.18f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(.1f, .1f, .7f));
    lightingShaderWithTexture.setMat4("model", model1 * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.1, 0.1, 0.44));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


}


void drawtrain(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO) {
    // Draw rail line

    glLineWidth(10.0f);
    
    drawCircleflat(2.7f, -.80f, -5.0f, 1.9f, 36);  // Larger circle behind
    drawCircleflat(2.7f, -.80f, -5.0f, 1.5f, 36);   // Larger circle in front

    for (int i = 0; i < 60; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 60;
        GLfloat x1 = 1.5f * cos(angleRad)+2.7, z1 = 1.5f * sin(angleRad)-5.0, y1 = -.80f;
        GLfloat x2 = 1.9f * cos(angleRad)+2.7, z2 = 1.9f * sin(angleRad)-5.0, y2 = -.80f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    static float trainAngle = 0.0f;

    // In your main loop:
    trainAngle += 0.5f; // Increase the train's angle by 0.5 degrees every frame
    if (trainAngle > 360.0f) {
        trainAngle -= 360.0f; // Reset the train's angle when it reaches 360 degrees
    }

    glm::mat4 trainModel = glm::mat4(1.0f);
    trainModel = glm::translate(trainModel, glm::vec3(1.7f * cos(glm::radians(trainAngle))+2.7, 0.0f+0.2, 1.7f * sin(glm::radians(trainAngle)))); // move to the correct position on the circular path
    trainModel = glm::scale(trainModel, glm::vec3(0.3f, 0.5f, 0.4f));
    lightingShaderWithTexture.setMat4("model", trainModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(0.4, 0.8, 0.3));
    glBindVertexArray(cubeVAO);
    //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
    glm::vec3 seatPosition = glm::vec3(1.7f * cos(glm::radians(trainAngle + 35.0f))+2.7, 0.12f+0.6, 1.7f * sin(glm::radians(trainAngle + 35.0f))-5.0f);
    //drawSeat(seatPosition, glm::vec3(0.5f, 0.5f, 0.5f), 0.0f, glm::vec3(1, 1, 0), glm::mat4(1.0f), glm::mat4(1.0f), lightingShaderWithTexture, cubeVAO);

}


void drawSeat(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO) {
    // Apply transformations
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 trainmodel = glm::mat4(1.0f);
    translation =glm::vec3 (2.0,0.0,-5.0);
    model = glm::translate(model, translation); // Apply translation
    // Apply rotation
    model = glm::scale(model, scale);
    model1 = glm::scale(model1, scale);

    model = model * alTogether;
    model1 = model1 * alTogether2;
    lightingShaderWithTexture.setMat4("model", model);

    lightingShaderWithTexture.setVec3("curveColor", color);


    // Set line width

    glLineWidth(10.0f);

    // Draw circles
    drawCircle(0.0f, -0.2f, 0.0f, 0.2f, 36);  // Small circle behind
    drawCircle(0.0f, 0.2f, 0.0f, 0.2f, 36);   // Large circle in front
    drawCircle(0.4f, -0.2f, 0.0f, 0.2f, 36);  // Larger circle behind
    drawCircle(0.4f, 0.2f, 0.0f, 0.2f, 36);   // Larger circle in front
    glm::mat4 seatModel = glm::mat4(1.0f);
    //seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
    //seatModel = glm::rotate(seatModel, glm::radians(-rotationAngle), glm::vec3(0.0f, 0.0f, 1.0f)); // rotate to match the Ferris wheel's rotation
    seatModel = glm::translate(seatModel, glm::vec3(-0.25f, 0.1f, -0.28f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(.95f, 0.5f, .5f));
    
    lightingShaderWithTexture.setMat4("model", model * seatModel );
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 0));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


}

void drawCone(Shader& lightingShader, glm::vec3 color, glm::vec3 trans, glm::vec3 scale, glm::vec3 rot, float radius, float height)
{
    lightingShader.use();
    // building model matrix
    glm::mat4 model(1.0f);
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix;
    translateMatrix = glm::translate(model, trans);
    rotateXMatrix = glm::rotate(translateMatrix, glm::radians(rot[0]), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(rotateXMatrix, glm::radians(rot[1]), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(rotateYMatrix, glm::radians(rot[2]), glm::vec3(0.0f, 0.0f, 1.0f));
    model = glm::scale(rotateZMatrix, scale);


    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;
    int n = 36;

    for (int i = 0; i < n; i++)
    {
        float angle = 2.0f * 3.1416 * i / n;
        float x1 = radius * cos(angle), y1 = radius * sin(angle), z1 = -height;
        vertices.push_back(x1); vertices.push_back(y1); vertices.push_back(z1);
        vertices.push_back(color.r); vertices.push_back(color.g); vertices.push_back(color.b);

    }
    vertices.push_back(0); vertices.push_back(0); vertices.push_back(height);

    for (int i = 0; i < n - 1; i++)
    {
        indices.push_back(i);
        indices.push_back(n);
        indices.push_back(i + 1);
        indices.push_back(i + 1);
        indices.push_back(n);
        indices.push_back(i + 2);
    }
    indices.push_back(n - 1);
    indices.push_back(n);
    indices.push_back(0);
    indices.push_back(0);
    indices.push_back(n);
    indices.push_back(1);

    unsigned int cubeVAO, cubeVBO, cubeEBO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);


    // setting up materialistic property
    lightingShader.setVec3("material.ambient", color);
    lightingShader.setVec3("material.diffuse", color);
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);

    lightingShader.setMat4("model", model);


    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

    glBindVertexArray(0);
    glDeleteVertexArrays(1, &cubeVAO);
    glDeleteBuffers(1, &cubeVBO);
    glDeleteBuffers(1, &cubeEBO);
}

void shop(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& shop1roof, Cube& shop1stick, Cube& shop1base, Cube& pastry)
{


    //roof
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 translate2,rotat;
    rotat = glm::mat4(1.0f);

    scale = glm::scale(model, glm::vec3(1.8, 0.3, -1.1));
    translate = glm::translate(model, glm::vec3(.0, 2.0, 2.8));
    float angle = glm::radians(90.0f);
    rotat = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, (float)255 / 255, (float)165 / 255, (float)0 / 255);
    shaderActivate(lightingShaderWithTexture);
    shop1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //lower stick 
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.1, -.1, -1.4));
    translate = glm::translate(model, glm::vec3(-.10, 8.0, 2.0));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, (float)139 / 255, (float)69 / 255, (float)19 / 255);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //Base
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.8, 0.7, -.6));
    translate = glm::translate(model, glm::vec3(.0, -1.2, 5.2));
    angle = glm::radians(90.0f);
    rotat= glm::rotate(glm::mat4(1.0f), angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, (float)255 / 255, (float)165 / 255, (float)0 / 255);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    //Base front
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.8, 0.7, -.01));
    translate = glm::translate(model, glm::vec3(.0, -1.2, 311.0));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, (float)255 / 255, (float)165 / 255, (float)0 / 255);
    pastry.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, -.9, .1));
    translate = glm::translate(model, glm::vec3(-3.3, .750, -.10));
    model = alTogether * translate * scale;
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, -.9, .1));
    translate = glm::translate(model, glm::vec3(-3.3, .750, -1.80));
    model = alTogether * translate * scale;
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //Left leg 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, -1.5, .1));
    translate = glm::translate(model, glm::vec3(-4.13, .750, -.10));
    model = alTogether * translate * scale;
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    

    //Right leg 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, -1.5, .1));
    translate = glm::translate(model, glm::vec3(-4.13, .750, -1.80));
    model = alTogether * translate * scale;
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);
}

void shop1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& shop1roof, Cube& shop1stick, Cube& shop1base)
{
    ////basefront
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotat = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, 1.45, 0.05));
    translate = glm::translate(model, glm::vec3(2.5, -.6750, -76.5));
    float angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.772, 0.741, 0.486);
    shaderActivate(lightingShaderWithTexture);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    ////base
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.350, 0.1, 1.0));
    translate = glm::translate(model, glm::vec3(1.76, -9.0, -3.83));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.537, 0.403, 0.470);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    ////baseback
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, .5, 0.05));
    translate = glm::translate(model, glm::vec3(2.5, -1.80, -60.0));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.811, 0.466, 0.643);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    ////baseleft
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.5, .84));
    translate = glm::translate(model, glm::vec3(50.0, -1.80, -4.50));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.811, 0.466, 0.643);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    ////baseright
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.5, .84));
    translate = glm::translate(model, glm::vec3(69.0, -1.80, -4.50));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.811, 0.466, 0.643);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    ////stick1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(25.50,-0.64,-31.5));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.788, 0.113, 0.458);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    ////stick2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(33.50, -0.64, -31.5));
    angle = glm::radians(90.0f);
    rotat = glm::rotate(model, angle, glm::vec3(.0f, 1.0f, 0.0f));
    model = alTogether * rotat * scale * translate;
    ////drawCube(cubeVAO, lightingShader, model, 0.788, 0.113, 0.458);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, 0.05, 1.2));
    translate = glm::translate(model, glm::vec3(-5.3, 56.2, -3.0));
    angle = glm::radians(40.0f); // Replace 45.0f with the desired angle in degrees
    glm::mat4 rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(.0f, 0.0f, 1.0f));
    model = alTogether * rotateX * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.984, 0.835, 0.913);
    shop1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, 0.05, 1.2));
    translate = glm::translate(model, glm::vec3(-6.36, -32.5, -3.0));
    angle = glm::radians(-40.0f); // Replace 45.0f with the desired angle in degrees
    rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(.0f, 0.0f, 1.0f));
    model = alTogether * rotateX * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.984, 0.835, 0.913);
    shop1roof.drawCubeWithTexture(lightingShaderWithTexture, model);
}

void shop2(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& shop1roof, Cube& shop1stick, Cube& shop1base)
{

    float baseHeight = 1.0f;
    float width = 1.3f;
    float length = 0.4f;

    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, -10.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.0, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.772, 0.741, 0.486);
    shaderActivate(lightingShaderWithTexture);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    //base
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.450, 0.1, 1.0));
    translate = glm::translate(model, glm::vec3(-0.5, -10.0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.537, 0.403, 0.470);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseback
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, .5, 0.05));
    translate = glm::translate(model, glm::vec3(-0.5, -10.0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.811, 0.466, 0.643);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseleft
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.7, 0.5 + .24));
    translate = glm::translate(model, glm::vec3(-0.5, -10.0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.811, 0.466, 0.643);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseright
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.7, 0.5 + .24));
    translate = glm::translate(model, glm::vec3(-0.5, -10.0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.811, 0.466, 0.643);
    shop1base.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(-0.5, -10.0, -0.5));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.788, 0.113, 0.458);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(2 + 2 + .0, 0, -2 + 1));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.788, 0.113, 0.458);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(2 + 2 + .0, 0, -2 - 4 - 1));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.788, 0.113, 0.458);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(-5.0, 0, -2 - 4 - 1));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.788, 0.113, 0.458);
    shop1stick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.2, 0.05, 1.1 / 2.0));
    translate = glm::translate(model, glm::vec3(-0.5, 20 + 1, -2.5));
    float angle = glm::radians(40.0f); // Replace 45.0f with the desired angle in degrees
    glm::mat4 rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * rotateX * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.984, 0.835, 0.913);
    shop1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.2, 0.05, 1.1 / 2.0));
    translate = glm::translate(model, glm::vec3(-0.5, 20 + 10, .5));
    angle = glm::radians(-40.0f); // Replace 45.0f with the desired angle in degrees
    rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * rotateX * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.984, 0.835, 0.913);
    shop1roof.drawCubeWithTexture(lightingShaderWithTexture, model);
}

void dfloor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& floor)
{
    float baseHeight = 0.04;
    float width = 11.0;
    float length = 14.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 1, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(-0.5, -23.80, -.80));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    floor.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
}

void f_dining_room_wall(unsigned int& cubeVAO, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube,glm::vec3 color)
{
    float baseHeight = .7;
    float width = 0.01;
    float length = 13.0;
    glm::mat4 identityMatrix = glm::mat4(1.0f);

    // right
    float baseHeight1 = baseHeight / 3;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(identityMatrix, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(identityMatrix, glm::vec3(5.3, 2 * baseHeight1-1.41f, -11.70));
    model = alTogether  * translate * scale*rotateYMatrix;
    shaderActivate(lightingShaderWithTexture);
    lightingShaderWithTexture.setVec3("curveColor",color);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(identityMatrix, glm::vec3(width + 0.1 , baseHeight, length + 0.1));
    translate = glm::translate(identityMatrix, glm::vec3(-5.50, 2 * baseHeight1 - 1.41f, -11.70));
    model = alTogether * rotateYMatrix * translate * scale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    //back
    model = glm::mat4(1.0f);
    scale = glm::scale(identityMatrix, glm::vec3(width + 0.1, baseHeight, length + 0.1-2.20));
    translate = glm::translate(identityMatrix, glm::vec3(11.60, 2 * baseHeight1 - 1.41f, -5.50));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateYMatrix * translate * scale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front1
    model = glm::mat4(1.0f);
    scale = glm::scale(identityMatrix, glm::vec3(width + 0.1, baseHeight, 5.7 + 0.1 - 2.20));
    translate = glm::translate(identityMatrix, glm::vec3(-1.5, 2 * baseHeight1 - 1.41f, -5.50));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateYMatrix * translate * scale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front2
    model = glm::mat4(1.0f);
    scale = glm::scale(identityMatrix, glm::vec3(width + 0.1, baseHeight, 5.5 + 0.1 - 2.20));
    translate = glm::translate(identityMatrix, glm::vec3(-1.5, 2 * baseHeight1 - 1.41f, 2.20-0.2));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateYMatrix * translate * scale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    
}

void gate(unsigned int& cubeVAO, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Cube& cube1, glm::vec3 color1, glm::vec3 color2)
{
    float baseHeight = 3.3;
    float width = 0.4;
    float length = .4;
    glm::mat4 identityMatrix = glm::mat4(1.0f);

    // left
    float baseHeight1 = baseHeight / 3;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateYMatrix = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(identityMatrix, glm::vec3(width, baseHeight, length ));
    translate = glm::translate(identityMatrix, glm::vec3(0.0-1.9,  baseHeight1 - 2.0f, 1.30));
    model = alTogether * translate * scale * rotateYMatrix;
    shaderActivate(lightingShaderWithTexture);
    lightingShaderWithTexture.setVec3("curveColor", color1);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right
    scale = glm::scale(identityMatrix, glm::vec3(width, baseHeight, length));
    translate = glm::translate(identityMatrix, glm::vec3(1.60, baseHeight1 - 2.0f, 1.30));
    model = alTogether * translate * scale * rotateYMatrix;
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    //banner
    width = 3.1;
    baseHeight = 1.0;
    length = .10;
    scale = glm::scale(identityMatrix, glm::vec3(width, baseHeight, length));
    translate = glm::translate(identityMatrix, glm::vec3(0.0 - 1.5, baseHeight1 , 1.30));
    model = alTogether * translate * scale * rotateYMatrix;
    shaderActivate(lightingShaderWithTexture);
    cube1.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    
}


void f_drawing_tv(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_tv, Cube& drawing_sound_box, Cube& drawing_cupboard)
{
    float baseHeight = 0.4;
    float width = 0.6;
    float length = 0.02;

    // tv
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(0.0, 0.5, 1.475));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    shaderActivate(lightingShaderWithTexture);
    drawing_tv.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // sound box
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 0.05, length));
    translate = glm::translate(model, glm::vec3(0.0, 0.4, 1.475));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    drawing_sound_box.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // cupboard
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 0.2, length + 0.15));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 1.475 - 0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    drawing_cupboard.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
}


void fan(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& b, Cube& c, float x, float y, float z)
{
    float bladel = 1.5;
    float bladew = 0.2;
    float bladeh = 0.01;

    // Center
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 scale2 = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.27, 0.3, 0.27));
    scale2 = glm::scale(model, glm::vec3(0.5, 0.5, 0.5));
    translate = glm::translate(model, glm::vec3(-0.67, 0.0, -0.4));
    translate2 = glm::translate(model, glm::vec3(0.0, 1.35, 0.0));
    translate3 = glm::translate(model, glm::vec3(x, y, z));
    model = alTogether * translate3 * translate2 * scale2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    c.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(45.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    b.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(165.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    b.drawCubeWithTexture(lightingShaderWithTexture, model);


    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(285.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    b.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz) {
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx, ty, tz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(sx, sy, sz));
    model = translateMatrix * scaleMatrix;
    return model;
}

void f_dining_jug(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, float tx, float ty, float tz, float sx, float sy, float sz)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // toilet
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(sx + 0.08, sy + 0.17, sz + 0.05));
    translate = glm::translate(model, glm::vec3(tx + 2.5 + 1.0, ty + 0.44, tz -2.1));
    translate2 = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model =  alTogether * translate * rotateM * translate2 * scale;
    //bucket->draw(lightingShaderWithTexture, model);
}

void getCurrentTime(int& hours, int& minutes, int& seconds) {
    time_t currentTime = time(nullptr); // Get current UNIX timestamp
    struct tm* timeinfo;
    timeinfo = localtime(&currentTime);
    seconds = timeinfo->tm_sec;
    minutes = timeinfo->tm_min;
    hours = timeinfo->tm_hour;
}
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.ProcessKeyboard(FORWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.ProcessKeyboard(LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.ProcessKeyboard(RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera.ProcessKeyboard(UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        camera.ProcessKeyboard(DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_RIGHT, deltaTime);
    }

    /*if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
    {
        if (rotateAxis_X) rotateAngle_X -= 0.1;
        else if (rotateAxis_Y) rotateAngle_Y -= 0.1;
        else rotateAngle_Z -= 0.1;
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) translate_Y += 0.001;*/
    //if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) translate_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) translate_X += 0.001;
    /*if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) translate_X -= 0.001;*/
    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) translate_Z += 0.01;
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) translate_Z -= 0.01;
    //if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) scale_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) scale_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) scale_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) scale_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) scale_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) scale_Z -= 0.001;

    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        rotateAngle_X += 10.0;
        rotateAxis_X = 1.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        rotateAngle_Y += 10.0;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 1.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    {
        rotateAngle_Z += 10.0;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 1.0;
    }

   
  


    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
        if(isCubeMovingUp) {
            isCubeMovingUp = false;
            isCubeMovingDown = false;
        }
        else
        {
            isCubeMovingUp = true;
            isCubeMovingDown = false;
        }
    }
    

    // fan
    static bool keyPressed = false;
    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS && !keyPressed)
    {
        fanSwitch = !fanSwitch; // Toggle the fanSwitch
        keyPressed = true;      // Mark the key as pressed
    }
    else if (glfwGetKey(window, GLFW_KEY_0) == GLFW_RELEASE)
    {
        keyPressed = false; // Reset the key state when released
    }


    // TV
    static bool shiftPressed = false;
    if (glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_PRESS && !shiftPressed)
    {
        channel += 1;
        channel %= 4;  // Ensure the channel cycles between 0 and 3
        shiftPressed = true;  // Mark the Shift key as pressed
    }
    else if (glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_RELEASE)
    {
        shiftPressed = false;  // Reset the key state when released
    }


    if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
    {
        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        drawing_light.turnAmbientOn();
        bed_room1_light.turnAmbientOn();
        dining_light.turnAmbientOn();
       
        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        drawing_light.turnDiffuseOff();
        bed_room1_light.turnDiffuseOff();
        dining_light.turnDiffuseOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        drawing_light.turnSpecularOff();
        bed_room1_light.turnSpecularOff();
        dining_light.turnSpecularOff();  
    }

    if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        drawing_light.turnDiffuseOn();
        bed_room1_light.turnDiffuseOn();
        dining_light.turnDiffuseOn();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        drawing_light.turnAmbientOff();
        bed_room1_light.turnAmbientOff();
        dining_light.turnAmbientOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        drawing_light.turnSpecularOff();
        bed_room1_light.turnSpecularOff();
        dining_light.turnSpecularOff();
        
    }

    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        drawing_light.turnDiffuseOff();
        bed_room1_light.turnDiffuseOff();
        dining_light.turnDiffuseOff();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        drawing_light.turnAmbientOff();
        bed_room1_light.turnAmbientOff();
        dining_light.turnAmbientOff();
        
        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        drawing_light.turnSpecularOn();
        bed_room1_light.turnSpecularOn();
        dining_light.turnSpecularOn();
    }

    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        drawing_light.turnDiffuseOn();
        bed_room1_light.turnDiffuseOn();
        dining_light.turnDiffuseOn();

        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        drawing_light.turnAmbientOn();
        bed_room1_light.turnAmbientOn();
        dining_light.turnAmbientOn();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        drawing_light.turnSpecularOn();
        bed_room1_light.turnSpecularOn();
        dining_light.turnSpecularOn();
    }

    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
    {
        //pointlight2.turnOn();
        drawing_light.turnOn();
        bed_room1_light.turnOn();
        dining_light.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
    {
        //pointlight2.turnOff();
        drawing_light.turnOff();
        bed_room1_light.turnOff();
        dining_light.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
    {
        pointlight3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    {
        pointlight3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    {
        pointlight1.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        pointlight1.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        if (isFlywheelRotating)
        {
            isFlywheelRotating = false;
        }
        else
        {
            isFlywheelRotating = true;
        }
    }

    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOn();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOn();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOn();
        //pointlight4.turnSpecularOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        /*cout << "1 " << pointlight1.isOn() << endl;
        cout << pointlight2.isOn() << endl;
        cout << pointlight3.isOn() << endl;*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOff();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOff();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOff();
        //pointlight4.turnSpecularOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}
void load_texture(unsigned int& texture, string image_name, GLenum format)
{
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    int width, height, nrChannels;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(image_name.c_str(), &width, &height, &nrChannels, 0);
    if (data)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    else
    {
        cout << "Failed to load texture " << image_name << endl;
    }
    stbi_image_free(data);
}

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureWrappingModeS);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureWrappingModeT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, textureFilteringModeMin);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, textureFilteringModeMax);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}
void stage(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Cube& stagewall)
{
    float baseHeight = 0.21;
    float width = 2.30;
    float length = .80;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(-.65, -0.9, -2.55));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    width = 3.0;
    length = 0.01;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.30, 1.8, .01));
    translate = glm::translate(model, glm::vec3(-1.50, -0.20, -4.60));
    model = alTogether * translate2 * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    stagewall.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


