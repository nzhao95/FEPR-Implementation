#version 450 core // Minimal GL version support expected from the GPU

struct LightSource {
    vec3 position;
    vec3 color;
    float intensity;
    int isActive;
    mat4 depthMVP;
};

int number_of_lights = 3;

uniform LightSource lightSources[ 3 ];

uniform sampler2D lightShadowMaps[ 3 ];


struct Material {
    vec3 albedo;
    float shininess;
    vec3 F0;
    float roughness;
    float metalness;
};

uniform Material material;
uniform sampler2D albedoTex;

uniform int model;


in vec3 fPosition; // Shader input, linearly interpolated by default from the previous stage (here the vertex shader)
in vec3 fPositionWorldSpace;
in vec3 fNormal;
in vec2 fTexCoord;

out vec4 colorResponse; // Shader output: the color response attached to this fragment


uniform mat4 projectionMat, modelViewMat, normalMat;

float pi = 3.1415927;

float shadowCalculation(vec4 fragPosLightSpace, int light) {
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;
    float closestDepth = texture(lightShadowMaps[light], projCoords.xy).r;
    float currentDepth = projCoords.z;
    float shadow = currentDepth - 0.001 > closestDepth ? 1.0 : 0.0;
    //if ((projCoords-0.5).x * (projCoords-0.5).x + (projCoords-0.5).y * (projCoords-0.5).y>0.05){
    //    shadow = 1.0;
    //}
    return shadow;
}

float distanceAttenuation(vec4 fragPosLightSpace, int light) {
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    float d = projCoords.z;
    float ac = 1.0;
    float al = 1.0;
    float aq = 1.0;
    return 1/(ac + al * d + aq * d * d);
}

vec3 fresnelSchlick(float cosTheta, vec3 F0 ) {
  return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

float distribution (vec3 n, vec3 m, float roughness){
  return (roughness + 2)/(2*pi)*pow(dot(n,m), roughness);
}

float cookTorrance (vec3 n, vec3 wh, vec3 wi, vec3 wo) {
  return min(1.0, min(2*dot(n,wh)*dot(n,wi)/dot(wo,wh), 2*dot(n,wh)*dot(n,wo)/dot(wo,wh)));
}

void main() {
    vec3 n = normalize(fNormal);

    // Linear barycentric interpolation does not preserve unit vectors
    vec3 wo = normalize (-fPosition); // unit vector pointing to the camera
    vec3 radiance = vec3(0,0,0);

    if( dot( n , wo ) >= 0.0 ) {
        {
            for (int light = 0; light < number_of_lights; light++) {
              if( lightSources[light].isActive == 1 ) { // WE ONLY CONSIDER LIGHTS THAT ARE SWITCHED ON
                   {
                        vec3 wi = normalize ( vec3((modelViewMat * vec4(lightSources[light].position,1)).xyz) - fPosition ); // unit vector pointing to the light source (change if you use several light sources!!!)
                        if( dot( wi , n ) >= 0.0 ) { // WE ONLY CONSIDER LIGHTS THAT ARE ON THE RIGHT HEMISPHERE (side of the tangent plane)
                            vec3 wh = normalize( wi + wo ); // half vector (if wi changes, wo should change as well)
                            vec3 Li = lightSources[light].color * lightSources[light].intensity;

                            float shadow = shadowCalculation (lightSources[light].depthMVP * vec4(fPositionWorldSpace, 1.0), light);
                            float attenuation = distanceAttenuation (lightSources[light].depthMVP * vec4(fPositionWorldSpace, 1.0), light);

                            float cosTheta = dot(wh, wo);
                            vec3 F0 = mix(material.F0, material.albedo, material.metalness);

                            vec3 fresnelTerm = fresnelSchlick(cosTheta, F0);
                            float D = distribution(n, wh, material.roughness);
                            float divTerm = 1/(4*dot(n,wi)*dot(n,wo));
                            float G = cookTorrance(n, wh, wi, wo);

                            if (model == 0) {
                              radiance = radiance +
                                Li // light color
                                // * texture(albedoTex, fTexCoord).rgb
                                *material.albedo
                                * ( max(dot(n,wi),0.0) + pow(max(dot(n,wh),0.0),material.shininess) )
                                * (1 - shadow)
                                * attenuation
                                ;
                            }
                        }
                    }
                }
            }
        }
    }

    colorResponse = vec4 (radiance, 1.0); // Building an RGBA value from an RGB one.
}
