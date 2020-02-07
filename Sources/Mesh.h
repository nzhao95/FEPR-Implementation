#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <map>
#include <memory>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Transform.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Mesh : public Transform {
public:
	virtual ~Mesh ();

	inline const std::vector<glm::vec3> & vertexPositions () const { return m_vertexPositions; }
	inline std::vector<glm::vec3> & vertexPositions () { return m_vertexPositions; }

	inline const std::vector<glm::vec3> & vertexNormals () const { return m_vertexNormals; }
	inline std::vector<glm::vec3> & vertexNormals () { return m_vertexNormals; }

	inline const std::vector<glm::vec2> & vertexTexCoords () const { return m_vertexTexCoords; }
        inline std::vector<glm::vec2> & vertexTexCoords () { return m_vertexTexCoords; }

	inline const std::vector<glm::uvec3> & triangleIndices () const { return m_triangleIndices; }
	inline std::vector<glm::uvec3> & triangleIndices () { return m_triangleIndices; }

	/// Compute the parameters of a sphere which bounds the mesh
	void computeBoundingSphere (glm::vec3 & center, float & radius) const;

        void recomputePerVertexNormals (bool angleBased = false);
        void recomputePerVertexTextureCoordinates ( );
				void updateVelocity(float dt);
				void updatePositions();
				void internalForces();
				glm::vec3 externalForces(glm::vec3 x, float t);
				glm::vec3 elasticForce(glm::vec3 x, glm::vec3 x0);
				glm::vec3 dxExternalForces(glm::vec3 x);
				Eigen::VectorXf constraint ();
				void FEPR (float t);
				Eigen::VectorXf computeConstraint (Eigen::VectorXf q, Eigen::VectorXf qn_1);
				glm::vec3 computeVelocity(glm::vec3 x, glm::vec3 x0, glm::vec3 v, float t) ;
				Eigen::MatrixXf jacobian_c (Eigen::VectorXf q, Eigen::VectorXf qn_1);
				void updateJacobian (Eigen::MatrixXf &j, Eigen::VectorXf q);

	void init ();
	void render ();
	void clear ();

        void addPlan(float square_half_side = 1.2f);

private:
	std::vector<glm::vec3> m_vertexPositions;
	std::vector<glm::vec3> m_vertexNormals;
	std::vector<glm::vec2> m_vertexTexCoords;
	std::vector<glm::uvec3> m_triangleIndices;
	std::vector<glm::vec3> m_vertexVelocities;
	std::vector<glm::vec3> m_vertexStretch;
	std::vector<glm::vec3> m_vertexdxStretch;
	std::vector<glm::vec3> m_vertexRestPositions;
	std::vector<float> m_vertexRestLength;
	std::vector<glm::vec3> m_forceField;

	glm::vec3 barycenter0 = glm::vec3(0, 0, 0);
	glm::vec3 barycenter = glm::vec3(0, 0, 0);
	glm::vec3 barycenter_velocity = glm::vec3(0, 0, 0);
	glm::vec3 dbv = glm::vec3(0, 0, 0);
	glm::vec3 m_O = glm::vec3(0, 0, 0);

	Eigen::DiagonalMatrix<float, Eigen::Dynamic> m_D;

	GLuint m_vao = 0;
	GLuint m_posVbo = 0;
	GLuint m_normalVbo = 0;
	GLuint m_texCoordVbo = 0;
	GLuint m_ibo = 0;

	float m_mass = 1.f;
	float m_strech = 1.f;
	float m_ground = 0;
	float m_f = 1.f;
	float m_h = 0.05;
	float m_k = 1.f;
	float s = 1.f;
	float t = 1.f;
	float epsilon = 0.001;
};

#endif // MESH_H
