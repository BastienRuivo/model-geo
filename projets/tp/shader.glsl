#version 330

#ifdef VERTEX_SHADER

layout(location= 0) in vec3 position;
layout(location= 1) in vec2 texcoord;
layout(location= 2) in vec3 normal;
layout(location= 4) in int material;

uniform mat4 mvpMatrix;
uniform mat4 modelMatrix;

out vec3 pos;
out vec3 norm;
flat out int matId;
out vec2 uvs;

void main( )
{
    gl_Position = mvpMatrix * vec4(position, 1);

    pos = (modelMatrix * vec4(position, 1)).xyz;
    norm = ((modelMatrix * vec4(normal, 1)).xyz);
    
    matId = material;
    uvs = texcoord;
}
#endif

#ifdef FRAGMENT_SHADER
// doit calculer la couleur du fragment
in vec3 pos;
in vec3 norm;
flat in int matId;
in vec2 uvs;

const int MAX_MATERIALS = 8;
uniform vec4[MAX_MATERIALS] materials;
uniform sampler2D base_texture;    // declare une texture 2d

uniform vec3 source;
uniform vec4 lightColor;
uniform vec3 camera;
uniform float k;// constante de modèle de reflexion
uniform float shininess;
uniform float lightPower;

const float PI= 3.14159265359; 

void main( )
{

    vec3 camDir = normalize(camera - pos);
    vec3 lightDir = (source - pos);
    vec3 halfDir = normalize(camDir + lightDir);
    vec3 nn = normalize(norm);

    float cos_0 = max(0.0, dot(nn, normalize(lightDir)));
    float cos_0_half = max(0.0, dot(nn, halfDir));

    float fr = (shininess + 8) / (8 * PI) * pow(cos_0_half, shininess);

    vec3 color = lightColor.rgb * materials[matId].rgb * fr * cos_0 * lightPower / length(lightDir);

    gl_FragColor = vec4(color, 1.0);
}
#endif
