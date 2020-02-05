#define _USE_MATH_DEFINES

#include "Mesh.h"

#include <cmath>
#include <algorithm>

using namespace std;

Mesh::~Mesh () {
	clear ();
}

void Mesh::computeBoundingSphere (glm::vec3 & center, float & radius) const {
	center = glm::vec3 (0.0);
	radius = 0.f;
	for (const auto & p : m_vertexPositions)
		center += p;
	center /= m_vertexPositions.size ();
	for (const auto & p : m_vertexPositions)
		radius = std::max (radius, distance (center, p));
}

void Mesh::recomputePerVertexNormals (bool angleBased) {
	m_vertexNormals.clear ();
	// Change the following code to compute a proper per-vertex normal
        m_vertexNormals.resize (m_vertexPositions.size (), glm::vec3 (0.0, 0.0, 0.0));

        for( unsigned int tIt = 0 ; tIt < m_triangleIndices.size() ; ++tIt ) {
            glm::uvec3 t = m_triangleIndices[ tIt ];
            glm::vec3 n_t = glm::cross( m_vertexPositions[t[1]] - m_vertexPositions[t[0]] , m_vertexPositions[t[2]] - m_vertexPositions[t[0]] );
            m_vertexNormals[ t[0] ] += n_t;
            m_vertexNormals[ t[1] ] += n_t;
            m_vertexNormals[ t[2] ] += n_t;
        }
        for( unsigned int nIt = 0 ; nIt < m_vertexNormals.size() ; ++nIt ) {
            glm::normalize( m_vertexNormals[nIt] );
        }
}

void Mesh::recomputePerVertexTextureCoordinates() {
    m_vertexTexCoords.clear ();
    // Change the following code to compute a proper per-vertex texture coordinates
    m_vertexTexCoords.resize (m_vertexPositions.size (), glm::vec2 (0.0, 0.0));

    float xMin = FLT_MAX , xMax = FLT_MIN;
    float yMin = FLT_MAX , yMax = FLT_MIN;
    for ( glm::vec3 & p : m_vertexPositions) {
        xMin = std::min( xMin , p[0] );
        xMax = std::max( xMax , p[0] );
        yMin = std::min( yMin , p[1] );
        yMax = std::max( yMax , p[1] );
    }
    for( unsigned int pIt = 0 ; pIt < m_vertexTexCoords.size() ; ++pIt ) {
        m_vertexTexCoords[ pIt ] = glm::vec2( (m_vertexPositions[pIt][0] - xMin)/(xMax-xMin) , (m_vertexPositions[pIt][1] - yMin)/(yMax-yMin) );
    }
}

void Mesh::addPlan(float square_half_side) {
    m_vertexPositions.push_back(glm::vec3 (-square_half_side,-square_half_side , -1));
    m_vertexPositions.push_back(glm::vec3 (+square_half_side,-square_half_side , -1));
    m_vertexPositions.push_back(glm::vec3 (+square_half_side,+square_half_side , -1));
    m_vertexPositions.push_back(glm::vec3 (-square_half_side,+square_half_side , -1));

    m_vertexTexCoords.push_back(glm::vec2 (0.0, 0.0));
    m_vertexTexCoords.push_back(glm::vec2 (1.0, 0.0));
    m_vertexTexCoords.push_back(glm::vec2 (1.0, 1.0));
    m_vertexTexCoords.push_back(glm::vec2 (0.0, 1.0));

    m_vertexNormals.push_back(glm::vec3 (0,0, 1));
    m_vertexNormals.push_back(glm::vec3 (0,0, 1));
    m_vertexNormals.push_back(glm::vec3 (0,0, 1));
    m_vertexNormals.push_back(glm::vec3 (0,0, 1));

    m_triangleIndices.push_back( glm::uvec3( m_vertexPositions.size() - 4 , m_vertexPositions.size() - 3 , m_vertexPositions.size() - 2 ) );
    m_triangleIndices.push_back( glm::uvec3( m_vertexPositions.size() - 4 , m_vertexPositions.size() - 2 , m_vertexPositions.size() - 1 ) );
}

void Mesh::init () {

	glCreateBuffers (1, &m_posVbo); // Generate a GPU buffer to store the positions of the vertices
	size_t vertexBufferSize = sizeof (glm::vec3) * m_vertexPositions.size (); // Gather the size of the buffer from the CPU-side vector
	glNamedBufferStorage (m_posVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data store on the GPU
	glNamedBufferSubData (m_posVbo, 0, vertexBufferSize, m_vertexPositions.data ()); // Fill the data store from a CPU array

	glCreateBuffers (1, &m_normalVbo); // Same for normal
	glNamedBufferStorage (m_normalVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData (m_normalVbo, 0, vertexBufferSize, m_vertexNormals.data ());

	glCreateBuffers (1, &m_texCoordVbo); // Same for texture coordinates
	size_t texCoordBufferSize = sizeof (glm::vec2) * m_vertexTexCoords.size ();
	glNamedBufferStorage (m_texCoordVbo, texCoordBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData (m_texCoordVbo, 0, texCoordBufferSize, m_vertexTexCoords.data ());

	glCreateBuffers (1, &m_ibo); // Same for the index buffer, that stores the list of indices of the triangles forming the mesh
	size_t indexBufferSize = sizeof (glm::uvec3) * m_triangleIndices.size ();
	glNamedBufferStorage (m_ibo, indexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData (m_ibo, 0, indexBufferSize, m_triangleIndices.data ());

	glCreateVertexArrays (1, &m_vao); // Create a single handle that joins together attributes (vertex positions, normals) and connectivity (triangles indices)
	glBindVertexArray (m_vao);

	glEnableVertexAttribArray (0);
	glBindBuffer (GL_ARRAY_BUFFER, m_posVbo);
	glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof (GLfloat), 0);

	glEnableVertexAttribArray (1);
	glBindBuffer (GL_ARRAY_BUFFER, m_normalVbo);
	glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof (GLfloat), 0);

	glEnableVertexAttribArray (2);
	glBindBuffer (GL_ARRAY_BUFFER, m_texCoordVbo);
	glVertexAttribPointer (2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof (GLfloat), 0);

	glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, m_ibo);
	glBindVertexArray (0); // Desactive the VAO just created. Will be activated at rendering time.


	m_vertexVelocities.resize(m_vertexPositions.size(), glm::vec3(0.0));
	m_vertexStretch.resize(m_triangleIndices.size(), glm::vec3(0.0));
	m_vertexRestLength.resize(m_triangleIndices.size(), 0.f);
	m_vertexRestPositions = m_vertexPositions;
	for( unsigned int tIt = 0 ; tIt < m_triangleIndices.size() ; ++tIt ) {
		glm::uvec3 t = m_triangleIndices[ tIt ];
		const float delta_x1 = glm::length(m_vertexPositions[t[1]] - m_vertexPositions[t[0]]);
		const float delta_x2 = glm::length(m_vertexPositions[t[2]] - m_vertexPositions[t[0]]);
		m_vertexRestLength[tIt] = delta_x1+delta_x2;
	}
}

void Mesh::render () {
	size_t vertexBufferSize = sizeof (glm::vec3) * m_vertexPositions.size ();
	glNamedBufferSubData (m_posVbo, 0, vertexBufferSize, m_vertexPositions.data ()); //Update PosVBO with new positions
	glBindVertexArray (m_vao); // Activate the VAO storing geometry data
        glDrawElements (GL_TRIANGLES, static_cast<GLsizei> (m_triangleIndices.size () * 3), GL_UNSIGNED_INT, 0);
        // Call for rendering: stream the current GPU geometry through the current GPU program
}

void Mesh::clear () {
	m_vertexPositions.clear ();
	m_vertexNormals.clear ();
	m_vertexTexCoords.clear ();
	m_triangleIndices.clear ();
	if (m_vao) {
		glDeleteVertexArrays (1, &m_vao);
		m_vao = 0;
	}
	if(m_posVbo) {
		glDeleteBuffers (1, &m_posVbo);
		m_posVbo = 0;
	}
	if (m_normalVbo) {
		glDeleteBuffers (1, &m_normalVbo);
		m_normalVbo = 0;
	}
	if (m_texCoordVbo) {
		glDeleteBuffers (1, &m_texCoordVbo);
		m_texCoordVbo = 0;
	}
	if (m_ibo) {
		glDeleteBuffers (1, &m_ibo);
		m_ibo = 0;
	}
}

glm::vec3 Mesh::externalForces(glm::vec3 x, float t) {
	float r = 0.2;
	glm::vec3 gravity = glm::vec3(0, -9.8, 0);
	glm::vec3 floor = glm::vec3(0, 9.8, 0);
	glm::vec3 force = glm::vec3(0.1f, 0, 0.6f + cos(t) * x[0] * 0.1f) * abs(x[1] - m_ground);
	if (t<0.2)
		return force;
	else return glm::vec3(0);
}

glm::vec3 Mesh::elasticForce(glm::vec3 x, glm::vec3 x0) { // consider the rhino as an elastic object
	return - m_k * m_f * (x - x0);
}


void Mesh::FEPR (float t) {
	const unsigned int n = m_vertexPositions.size();

	Eigen::VectorXf q (n * 6 + 2) ;
	Eigen::VectorXf qn_1 (n * 6 + 2) ;
	Eigen::VectorXf D_diagonal (n * 6 + 2);



	for (unsigned int nIt ; nIt < n; ++nIt) {
		glm::vec3 v = computeVelocity(nIt, t);
		glm::vec3 x = m_vertexPositions[nIt] + v * m_h;
		q(nIt*3) = x[0];
		q(nIt*3+1) = x[1];
		q(nIt*3+2) = x[2];
		q((nIt+n)*3) = v[0];
		q((nIt+n)*3+1) = v[1];
		q((nIt+n)*3+2) = v[2];
		D_diagonal(nIt*3) = 1; D_diagonal(nIt*3+1) = 1; D_diagonal(nIt*3+2) = 1;
		D_diagonal((nIt+n)*3) = pow(m_h,2); D_diagonal((nIt+n)*3+1) = pow(m_h,2); D_diagonal((nIt+n)*3+2) = pow(m_h,2);
	}

	D_diagonal(6*n) = 0.001; D_diagonal(6*n+1) = 0.001;

	Eigen::DiagonalMatrix<float, Eigen::Dynamic> D;
	D.diagonal() = D_diagonal;
	qn_1 = q;
	Eigen::VectorXf c = computeConstraint(q,  qn_1);

	unsigned int iterations = 0;

	while (iterations < 10 && c.sum() > 10e-7) {
		Eigen::MatrixXf jacobC = jacobian_c(q, qn_1);
		c = computeConstraint(q,  qn_1);
		Eigen::MatrixXf A = jacobC.transpose() * D * jacobC;

		Eigen::VectorXf lambda = A.colPivHouseholderQr().solve(-c);

		q = q - D.inverse() * jacobC * lambda ;

		iterations++;
	}

	for (unsigned int nIt ; nIt < n; ++nIt) {
		m_vertexPositions[nIt] = glm::vec3(q(nIt*3),q(nIt*3+1),q(nIt*3+2));
		m_vertexVelocities[nIt] = glm::vec3(q((nIt+n)*3),q((nIt+n)*3+1),q((nIt+n)*3+2));
	}

}

Eigen::VectorXf Mesh::computeConstraint (Eigen::VectorXf q, Eigen::VectorXf qn_1){
	Eigen::VectorXf c (7);
	float energy = 0.f;
	glm::vec3 linear_momentum (0);
	glm::vec3 angular_momentum (0);

	for (unsigned int nIt ; nIt < m_vertexPositions.size(); ++nIt) {
		unsigned int qIt = (nIt + m_vertexPositions.size()) * 3;
		const glm::vec3 x = glm::vec3(q(nIt*3), q(nIt*3+1), q(nIt*3+2));
		const glm::vec3 x0 = m_vertexRestPositions[nIt];
		const glm::vec3 v = glm::vec3(q(qIt), q(qIt+1), q(qIt+2));

		const glm::vec3 xn = m_vertexPositions[nIt];
		const glm::vec3 vn = m_vertexVelocities[nIt];

		const glm::vec3 xn_1 = glm::vec3(qn_1(nIt*3), qn_1(nIt*3+1), qn_1(nIt*3+2));
		const glm::vec3 vn_1 = glm::vec3(qn_1(qIt), qn_1(qIt+1), qn_1(qIt+2));

		energy += (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + pow(glm::length(x - x0), 2) * m_k)
		-(vn[0]*vn[0] + vn[1]*vn[1] + vn[2]*vn[2] + pow(glm::length(xn - x0), 2) * m_k); //Sum of cinetic and elastic energy
		linear_momentum += v - vn_1 + q(m_vertexPositions.size()*6) * (vn_1 - vn);
		angular_momentum += glm::cross(x,v) - glm::cross(xn_1, vn_1)
		+ q(m_vertexPositions.size()*6+1) * (glm::cross(xn_1, vn_1) - glm::cross(xn,vn));
	}

	c << energy, linear_momentum[0], linear_momentum[1], linear_momentum[2],
	angular_momentum[0], angular_momentum[1], angular_momentum[2];

	return c;
}

Eigen::MatrixXf Mesh::jacobian_c (Eigen::VectorXf q, Eigen::VectorXf qn_1) {
	unsigned int n = m_vertexPositions.size();
	Eigen::MatrixXf dc (n * 6 + 2, 7);
	glm::vec3 pn(0);
	glm::vec3 pn_1(0);
	glm::vec3 ln(0);
	glm::vec3 ln_1(0);
	for (unsigned int nIt ; nIt < n; ++nIt) {
		const glm::vec3 x = glm::vec3(q(nIt*3), q(nIt*3+1), q(nIt*3+2));
		const glm::vec3 x0 = m_vertexRestPositions[nIt];
		const glm::vec3 v = glm::vec3(q((nIt+n)*3), q((nIt+n)*3+1), q((nIt+n)*3+2));

		const glm::vec3 xn = m_vertexPositions[nIt];
		const glm::vec3 vn = m_vertexVelocities[nIt];

		const glm::vec3 xn_1 = glm::vec3(qn_1(nIt*3), qn_1(nIt*3+1), qn_1(nIt*3+2));
		const glm::vec3 vn_1 = glm::vec3(qn_1((nIt+n)*3), qn_1((nIt+n)*3+1), qn_1((nIt+n)*3+2));

		pn += vn;
		pn_1 += vn_1;
		ln += glm::cross(xn, vn);
		ln_1 += glm::cross(xn_1, vn_1);

		//Computing partial derivates according to xi, vi, s,t with 1 column for each condition

		dc(nIt*3, 0) = 2*(q(nIt*3)-x0[0]); dc(nIt*3+1, 0) = 2*(q(nIt*3+1)-x0[1]); dc(nIt*3+1, 0) = 2*(q(nIt*3+2)-x0[2]);
		dc((nIt+n)*3, 0) = 2*q((nIt+n)*3); dc((nIt+n)*3 + 1, 0) = 2*q((nIt+n)*3+1); dc((nIt+n)*3 + 2, 0) = 2*q((nIt+n)*3+2);

		dc(nIt*3, 1) = 0; dc(nIt*3+1, 1) = 0; dc(nIt*3+1, 1) = 0;
		dc(nIt*3, 2) = 0; dc(nIt*3+1, 2) = 0; dc(nIt*3+1, 2) = 0;
		dc(nIt*3, 3) = 0; dc(nIt*3+1, 3) = 0; dc(nIt*3+1, 3) = 0;

		dc((nIt+n)*3, 1) = 1; dc((nIt+n)*3 + 1, 1) = 0; dc((nIt+n)*3 + 2, 1) = 0;
		dc((nIt+n)*3, 2) = 0; dc((nIt+n)*3 + 1, 2) = 1; dc((nIt+n)*3 + 2, 2) = 0;
		dc((nIt+n)*3, 3) = 0; dc((nIt+n)*3 + 1, 3) = 0; dc((nIt+n)*3 + 2, 3) = 1;

		dc(nIt*3, 4) = 0; dc(nIt*3+1, 4) = v[2]; dc(nIt*3+1, 4) = -v[1];
		dc(nIt*3, 5) = -v[2]; dc(nIt*3+1, 5) = 0; dc(nIt*3+1, 5) = v[0];
		dc(nIt*3, 6) = v[1]; dc(nIt*3+1, 6) = -v[0]; dc(nIt*3+1, 6) = 0;

		dc((nIt+n)*3, 4) = 0; dc((nIt+n)*3 + 1, 4) = -x[2]; dc((nIt+n)*3 + 2, 4) = x[1];
		dc((nIt+n)*3, 5) = x[2]; dc((nIt+n)*3 + 1, 5) = 0; dc((nIt+n)*3 + 2, 5) = -x[0];
		dc((nIt+n)*3, 6) = -x[1]; dc((nIt+n)*3 + 1, 6) = x[0]; dc((nIt+n)*3 + 2, 6) = 0;

	}

	glm::vec3 ds = pn_1 - pn;
	glm::vec3 dt = ln_1 - ln;


	dc(6*n, 0) = 0; dc(6*n+1, 0) = 0;
	dc(6*n, 1) = ds[0]; dc(6*n+1, 1) = 0;
	dc(6*n, 2) = ds[1]; dc(6*n+1, 2) = 0;
	dc(6*n, 3) = ds[2]; dc(6*n+1, 3) = 0;
	dc(6*n, 4) = 0; dc(6*n+1, 4) = dt[0];
	dc(6*n, 5) = 0; dc(6*n+1, 5) = dt[1];
	dc(6*n, 6) = 0; dc(6*n+1, 6) = dt[2];
}


void Mesh::updateVelocity(float t) {
	for( unsigned int nIt = 0 ; nIt < m_vertexPositions.size() ; ++nIt ) {
		glm::vec3 x = m_vertexPositions[nIt];
		// m_vertexVelocities[nIt] = m_vertexVelocities[nIt] + m_vertexAccelerations[nIt] * m_h;
		m_vertexVelocities[nIt] += m_h * (externalForces(x,t) + elasticForce(x, m_vertexRestPositions[nIt])
		 + m_h * (- glm::vec3(m_k))
		 * m_vertexVelocities[nIt])/(glm::vec3(1) - (float) pow(m_h, 2) * glm::vec3(m_k));
		glm::vec3 angular_momentum = glm::cross(x - m_O, m_vertexVelocities[nIt]); // ||v||sin(theta).T
		glm::vec3 angular_vector = glm::normalize(glm::cross(x - m_O, angular_momentum));
		m_vertexVelocities[nIt] = glm::dot(m_vertexVelocities[nIt], angular_vector) * angular_vector;

	}
}

glm::vec3 Mesh::computeVelocity(unsigned int nIt, float t) {
	glm::vec3 x = m_vertexPositions[nIt];
	// m_vertexVelocities[nIt] = m_vertexVelocities[nIt] + m_vertexAccelerations[nIt] * m_h;
	glm::vec3 v = m_vertexVelocities[nIt] + m_h * (externalForces(x,t) + elasticForce(x, m_vertexRestPositions[nIt])
	 + m_h * (- glm::vec3(m_k))
	 * m_vertexVelocities[nIt])/(glm::vec3(1) - (float) pow(m_h, 2) * glm::vec3(m_k));
	glm::vec3 angular_momentum = glm::cross(x - m_O, m_vertexVelocities[nIt]); // ||v||sin(theta).T
	glm::vec3 angular_vector = glm::normalize(glm::cross(x - m_O, angular_momentum));
	v = glm::dot(v, angular_vector) * angular_vector;

	return v;
}

void Mesh::updatePositions() {

	for( unsigned int nIt = 0 ; nIt < m_vertexPositions.size() ; ++nIt ) {
		m_vertexPositions[nIt] = m_vertexPositions[nIt] + m_vertexVelocities[nIt] * m_h;
	}
	recomputePerVertexNormals();
}


glm::vec3 Mesh::computePosition(unsigned int nIt) {
		return m_vertexPositions[nIt] + m_vertexVelocities[nIt] * m_h;
}
