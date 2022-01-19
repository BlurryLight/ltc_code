#ifndef _LTC_
#define _LTC_


#include <glm/glm.hpp>
using namespace glm;

#include <iostream>
using namespace std;

struct LTC {

	// lobe magnitude
	float magnitude;

	// Average Schlick Fresnel term
	float fresnel;

	// parametric representation
	float m11, m22, m13;
	vec3 X, Y, Z;

	// matrix representation
	mat3 M;
	mat3 invM;
	float detM;
        float detInvM;

	LTC()
	{
		magnitude = 1;
		fresnel = 1;
		m11 = 1;
		m22 = 1;
		m13 = 0;
		X = vec3(1, 0, 0);
		Y = vec3(0, 1, 0);
		Z = vec3(0, 0, 1);
		update();
	}

	void update() // compute matrix from parameters
	{
          //注意glm是col major的
		M = mat3(X, Y, Z) *
			mat3(m11, 0, 0,
				0, m22, 0,
				m13, 0, 1);
		invM = inverse(M);
		detM = abs(glm::determinant(M));
                detInvM = abs(glm::determinant(invM));
	}

	float eval(const vec3& L) const
	{
                vec3 L_cosine = invM * L;
                vec3 Loriginal = normalize(L_cosine);//回到cosine分布的方向
                vec3 L_ = M * Loriginal;//再乘以m转到GGX方向，此时L_应该和L方向相同，只是有长度

		float l = length(L_);
                float l2 = length(L_cosine);
		float Jacobian = detM / (l*l*l);
                float Jacobian2 = detInvM / (l2 * l2 * l2);

		float D = 1.0f / 3.14159f * glm::max<float>(0.0f, Loriginal.z); //clamped cosine sample
		
//		float res = magnitude * D / Jacobian;
                float res = magnitude * D * Jacobian2;
		return res;
	}

        //看起来似乎是cosine weigted sampling
        // 不是
        // 返回的是M * L
        // 返回的是GGX分布的L

	vec3 sample(const float U1, const float U2) const
	{
		const float theta = acosf(sqrtf(U1));
		const float phi = 2.0f*3.14159f * U2;
		const vec3 L = normalize(M * vec3(sinf(theta)*cosf(phi), sinf(theta)*sinf(phi), cosf(theta)));
		return L;
	}
};

#endif
