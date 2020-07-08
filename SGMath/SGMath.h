/* 
	Simple Game Math

	Copyright (c) 2020 Patrick T. Coakley

	This software is provided 'as-is', without any express or implied
	warranty.  In no event will the authors be held liable for any damages
	arising from the use of this software.

	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:

	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution. 
*/

#ifndef SGMATH_H
#define SGMATH_H

#include <cmath>

// Forward Declaration
struct Matrix;
struct Quaternion;
struct Plane;

namespace sg
{
	// General math
	constexpr float Pi = 3.141593f;
	constexpr float TwoPi = 6.283185f;
	constexpr float HalfPi = 1.570796f;
	constexpr float Epsilon = 0.0000010f;
	constexpr float RadianToDegrees = 57.295780f;
	constexpr float DegreesToRadians = 0.017453f;

	constexpr float Lerp(float start, float end, float time) noexcept
	{
		return start + ((end - start) * time);
	}

	constexpr float Clamp(float min, float max, float f) noexcept
	{
		if (f < min)
		{
			return min;
		}

		if (f > max)
		{
			return max;
		}

		return f;
	}

	constexpr float ToRadians(float degrees) noexcept
	{
		return degrees * DegreesToRadians;
	}

	constexpr float ToDegrees(float radians) noexcept
	{
		return radians * RadianToDegrees;
	}

	inline float Sqrt(float f) noexcept
	{
		return sqrtf(f);
	}

	constexpr float Max(float a, float b) noexcept
	{
		if (a > b)
		{
			return a;
		}

		return b;
	}

	constexpr float Min(float a, float b) noexcept
	{
		if (a < b)
		{
			return a;
		}

		return b;
	}

	inline float Sin(float f) noexcept
	{
		return sinf(f);
	}

	inline float Cos(float f) noexcept
	{
		return cosf(f);
	}

	inline float Tan(float f) noexcept
	{
		return tanf(f);
	}

	inline float Asin(float f) noexcept
	{
		return asinf(f);
	}

	inline float Acos(float f) noexcept
	{
		return acosf(f);
	}

	inline float Atan(float f) noexcept
	{
		return atanf(f);
	}

	inline float Atan2(float a, float b) noexcept
	{
		return atan2f(a, b);
	}

	inline float Abs(float f) noexcept
	{
		return fabsf(f);
	}

	inline bool IsNaN(float f) noexcept
	{
		return isnan(f);
	}

	inline float Pow(float a, float b) noexcept
	{
		return powf(a, b);
	}

	inline float Random() noexcept
	{
		return static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
	}

	inline float Random(float max) noexcept
	{
		return max * Random();
	}

	inline float RandomRange(float min, float max) noexcept
	{
		return min + (max - min) * Random();
	}

	// 2D Vector definition
	struct Vector2
	{
		float X{};
		float Y{};

		Vector2() = default;
		Vector2(float x, float y);

		float& operator[](int i);
		const float& operator[](int i) const;
		Vector2 operator-() const;
		Vector2 operator+(const Vector2& vector2) const;
		Vector2 operator-(const Vector2& vector2) const;
		Vector2 operator*(const Vector2& vector2) const;
		Vector2 operator/(const Vector2& vector2) const;
		Vector2 operator*(float f) const;
		Vector2 operator/(float f) const;
		Vector2& operator+=(const Vector2& vector2);
		Vector2& operator-=(const Vector2& vector2);
		Vector2& operator+=(float f);
		Vector2& operator-=(float f);
		Vector2& operator*=(float f);
		Vector2& operator/=(float f);
		bool operator==(const Vector2& vector2) const;
		bool operator!=(const Vector2& vector2) const;
	};

	float Distance(const Vector2& lhs, const Vector2& rhs);
	float Dot(const Vector2& lhs, const Vector2& rhs);
	Vector2 Lerp(const Vector2& lhs, const Vector2& rhs, float t);
	float Magnitude(const Vector2& vector2);
	float MagnitudeSquared(const Vector2& vector2);
	Vector2 Normalized(const Vector2& vector2);
	Vector2 Reflect(const Vector2& vector2, const Vector2& normal);
	Vector2 Rotate(const Vector2& vector2, float r);
	Vector2 Scale(const Vector2& vector2, float f);
	Vector2 Transform(const Vector2& vector2, const Matrix& matrix);
	Vector2 operator*(float f, const Vector2& vector2);
	Vector2 operator/(float f, const Vector2& vector2);

	// 3D Vector definition
	struct Vector3
	{
		float X{};
		float Y{};
		float Z{};

		Vector3() = default;
		Vector3(float x, float y, float z);

		float& operator[](int i);
		const float& operator[](int i) const;
		Vector3 operator-() const;
		Vector3 operator+(const Vector3& vector3) const;
		Vector3 operator-(const Vector3& vector3) const;
		Vector3 operator*(const Vector3& vector3) const;
		Vector3 operator/(const Vector3& vector3) const;
		Vector3 operator*(float f) const;
		Vector3 operator/(float f) const;
		Vector3& operator+=(const Vector3& vector3);
		Vector3& operator-=(const Vector3& vector3);
		Vector3& operator+=(float f);
		Vector3& operator-=(float f);
		Vector3& operator*=(float f);
		Vector3& operator/=(float f);
		bool operator==(const Vector3& vector3) const;
		bool operator!=(const Vector3& vector3) const;
	};

	Vector3 Cross(const Vector3& lhs, const Vector3& rhs);
	float Distance(const Vector3& lhs, const Vector3& rhs);
	float Dot(const Vector3& lhs, const Vector3& rhs);
	Vector3 Lerp(const Vector3& lhs, const Vector3& rhs, float t);
	float Magnitude(const Vector3& vector3);
	float MagnitudeSquared(const Vector3& vector3);
	Vector3 Normalized(const Vector3& vector3);
	Vector3 Reflect(const Vector3& vector3, const Vector3& normal);
	Vector3 Rotate(const Vector3& vector3, const Quaternion& quaternion);
	Vector3 Scale(const Vector3& vector3, float f);
	Vector3 Transform(const Vector3& vector3, const Matrix& matrix);
	Vector3 operator*(float f, const Vector3& vector3);
	Vector3 operator/(float f, const Vector3& vector3);

	// 4D Vector definition
	struct Vector4
	{
		float X{};
		float Y{};
		float Z{};
		float W{};

		Vector4() = default;
		explicit Vector4(const Vector3& vector3);
		Vector4(float x, float y, float z, float w);

		float& operator[](int i);
		const float& operator[](int i) const;
		Vector4 operator-() const;
		Vector4 operator+(const Vector4& vector4) const;
		Vector4 operator-(const Vector4& vector4) const;
		Vector4 operator*(const Vector4& vector4) const;
		Vector4 operator/(const Vector4& vector4) const;
		Vector4 operator*(float f) const;
		Vector4 operator/(float f) const;
		Vector4& operator+=(const Vector4& vector4);
		Vector4& operator-=(const Vector4& vector4);
		Vector4& operator+=(float f);
		Vector4& operator-=(float f);
		Vector4& operator*=(float f);
		Vector4& operator/=(float f);
		bool operator==(const Vector4& vector4) const;
		bool operator!=(const Vector4& vector4) const;
	};

	float Distance(const Vector4& lhs, const Vector4& rhs);
	float Dot(const Vector4& lhs, const Vector4& rhs);
	Vector4 Lerp(const Vector4& lhs, const Vector4& rhs, float t);
	float Magnitude(const Vector4& vector4);
	float MagnitudeSquared(const Vector4& vector4);
	Vector4 Normalized(const Vector4& vector4);
	Vector4 Reflect(const Vector4& vector4, const Vector4& normal);
	Vector4 operator*(float f, const Vector4& vector4);
	Vector4 operator/(float f, const Vector4& vector4);

	// 4x4 Matrix definition
	struct Matrix
	{
		Matrix() = default;
		explicit Matrix(float f);
		explicit Matrix(float m[]);
		Matrix(float m00, float m01, float m02, float m03, float m04, float m05,
			float m06, float m07, float m08, float m09, float m10, float m11,
			float m12, float m13, float m14, float m15);

		float& operator[](unsigned int m);
		const float& operator[](unsigned int m) const;
		float& operator()(unsigned int col, unsigned int row);
		const float& operator()(unsigned int col, unsigned int row) const;

		Matrix operator*(float f) const;
		Matrix operator*(const Matrix& matrix) const;
		static const Matrix Identity;

	private:
		float mMatrix[16]{};
	};

	float Determinant(const Matrix& matrix);
	Matrix Inverse(const Matrix& matrix);
	Matrix Lerp(const Matrix& lhs, const Matrix& rhs, float t);
	Matrix LookAt(const Vector3& position, const Vector3& target,
		const Vector3& up);
	Matrix MakeOrthographic(float l, float r, float t, float b, float n, float f);
	Matrix MakeFrustum(float l, float r, float b, float t, float n, float f);
	Matrix MakePerspective(float fovy, float a, float n, float f);
	Matrix MakeRotationFromAxisAngle(float x, float y, float z, float r);
	Matrix MakeRotationFromAxisAngle(const Vector3& axis, float r);
	Matrix MakeRotationFromAxisAngle(const Matrix& matrix, const Vector3& axis,
		float r);
	Matrix MakeRotationFromEuler(float roll, float pitch, float yaw);
	Matrix MakeRotationFromEuler(const Vector3& euler);
	Matrix MakeRotationX(float r);
	Matrix MakeRotationY(float r);
	Matrix MakeRotationZ(float r);
	Matrix MakeScale(float x, float y, float z);
	Matrix MakeScale(const Vector3& vector3);
	Matrix MakeShear(float x, float y, float z);
	Matrix MakeShear(const Vector3& vector3);
	Matrix MakeTranslationMatrix(float x, float y, float z);
	Matrix MakeTranslationMatrix(const Vector3& vector3);
	Matrix Transposed(const Matrix& matrix);
	Matrix operator*(float f, const Matrix& matrix4);

	// Quaternion definition
	struct Quaternion
	{
		float X{};
		float Y{};
		float Z{};
		float W{};

		Quaternion() = default;
		Quaternion(float x, float y, float z, float w);
		Quaternion(const Vector3& vector3, float w);

		float& operator[](int i);
		const float& operator[](int i) const;
		Quaternion operator-() const;
		Quaternion operator+(const Quaternion& quaternion) const;
		Quaternion operator-(const Quaternion& quaternion) const;
		Quaternion operator*(const Quaternion& quaternion) const;
		Quaternion operator+(float f) const;
		Quaternion operator-(float f) const;
		Quaternion operator*(float f) const;
		Quaternion operator/(float f) const;
		Quaternion& operator+=(float f);
		Quaternion& operator-=(float f);
		Quaternion& operator*=(float f);
		Quaternion& operator/=(float f);
		bool operator==(const Quaternion& quaternion) const;
		bool operator!=(const Quaternion& quaternion) const;
		static const Quaternion Identity;
	};

	float AngleBetween(const Quaternion& lhs, const Quaternion& rhs);
	float Dot(const Quaternion& lhs, const Quaternion& rhs);
	Quaternion Lerp(const Quaternion& lhs, const Quaternion& rhs, float t);
	float Magnitude(const Quaternion& quaternion);
	float MagnitudeSquared(const Quaternion& quaternion);
	Quaternion NLerp(const Quaternion& lhs, const Quaternion& rhs, float t);
	Quaternion Normalized(const Quaternion& quaternion);
	Quaternion RotateTowards(const Quaternion& from, const Quaternion& to,
		float max);
	Quaternion MakeQuaternionFromAxisAngle(float x, float y, float z, float angle);
	Quaternion MakeQuaternionFromAxisAngle(const Vector3& axis, float angle);
	Quaternion MakeQuaternionFromEuler(float yaw, float pitch, float roll);
	Quaternion MakeQuaternionFromEuler(const Vector3& euler);
	Quaternion MakeQuaternionFromMatrix(const Matrix& matrix);
	Quaternion Slerp(const Quaternion& lhs, const Quaternion& rhs, float t);
	Matrix ToMatrix(const Quaternion& quaternion);
	Quaternion operator+(float f, const Quaternion& quaternion);
	Quaternion operator-(float f, const Quaternion& quaternion);
	Quaternion operator*(float f, const Quaternion& quaternion);
	Quaternion operator/(float f, const Quaternion& quaternion);

	// Ray definition
	struct Ray
	{
		Vector3 Origin{};
		Vector3 Direction{};

		Ray() = default;
		Ray(const Vector3& origin, const Vector3& direction)
			: Origin(origin), Direction(direction) {}

		Vector3 At(float t) const;
	};

	float Intersects(const Ray& ray, const Plane& plane);

	// Plane definition
	struct Plane
	{
		Vector3 Normal{};
		float Offset{};

		Plane() = default;
		Plane(const Vector3& normal, float offset)
			: Normal(normal), Offset(offset) {}
	};

} // namespace sg

#endif // SGMATH_H
