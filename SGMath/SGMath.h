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

#include <cassert>
#include <cmath>
#include <cstring>

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
		return a > b ? a : b;
	}

	constexpr float Min(float a, float b) noexcept
	{
		return a < b ? a : b;
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

		constexpr Vector2() noexcept = default;
		constexpr Vector2(float X, float Y) noexcept
			: X(X), Y(Y)
		{
		}

		float& operator[](int i) noexcept
		{
			assert(i >= 0 && i < 2 && "Out of range: Vector2::Operator[]\n");
			return (&X)[i];
		}

		constexpr const float& operator[](int i) const noexcept
		{
			assert(i >= 0 && i < 2 && "Out of range: Vector2::Operator[]\n");
			return (&X)[i];
		}

		constexpr Vector2 operator-() const noexcept
		{
			return Vector2(-X, -Y);
		}

		constexpr Vector2 operator+(const Vector2& vector2) const noexcept
		{
			return Vector2(X + vector2.X, Y + vector2.Y);
		}

		constexpr Vector2 operator-(const Vector2& vector2) const noexcept
		{
			return Vector2(X - vector2.X, Y - vector2.Y);
		}

		constexpr Vector2 operator*(const Vector2& vector2) const noexcept
		{
			float x = X * vector2.X;
			float y = Y * vector2.Y;

			return Vector2(x, y);
		}

		constexpr Vector2 operator/(const Vector2& vector2) const noexcept
		{
			float x = X / vector2.X;
			float y = Y / vector2.Y;

			return Vector2(x, y);
		}

		constexpr Vector2 operator*(float f) const noexcept
		{
			return Vector2(X * f, Y * f);
		}

		constexpr Vector2 operator/(float f) const noexcept
		{
			f = 1.0f / f;
			return Vector2(X * f, Y * f);
		}

		Vector2& operator+=(const Vector2& vector2) noexcept
		{
			X += vector2.X;
			Y += vector2.Y;

			return *this;
		}

		Vector2& operator-=(const Vector2& vector2) noexcept
		{
			X -= vector2.X;
			Y -= vector2.Y;

			return *this;
		}

		Vector2& operator+=(float f) noexcept
		{
			X += f;
			Y += f;

			return *this;
		}

		Vector2& operator-=(float f) noexcept
		{
			X -= f;
			Y -= f;

			return *this;
		}

		Vector2& operator*=(float f) noexcept
		{
			X *= f;
			Y *= f;

			return *this;
		}

		Vector2& operator/=(float f) noexcept
		{
			f = 1.0f / f;
			X *= f;
			Y *= f;

			return *this;
		}

		constexpr bool operator==(const Vector2& vector2) const noexcept
		{
			return X == vector2.X && Y == vector2.Y;
		}

		constexpr bool operator!=(const Vector2& vector2) const noexcept
		{
			return !(vector2 == *this);
		}
	};

	// 3D Vector definition
	struct Vector3
	{
		float X{};
		float Y{};
		float Z{};

		constexpr Vector3() noexcept = default;
		constexpr Vector3(float x, float y, float z) noexcept
			: X(x), Y(y), Z(z)
		{
		}

		float& operator[](int i) noexcept
		{
			assert(i >= 0 && i < 3 && "Out of range: Vector3::Operator[]\n");
			return (&X)[i];
		}

		constexpr const float& operator[](int i) const noexcept
		{
			assert(i >= 0 && i < 3 && "Out of range: Vector3::Operator[]\n");
			return (&X)[i];
		}

		constexpr Vector3 operator-() const noexcept
		{
			return Vector3(-X, -Y, -Z);
		}

		constexpr Vector3 operator+(const Vector3& vector3) const noexcept
		{
			return Vector3(X + vector3.X, Y + vector3.Y, Z + vector3.Z);
		}

		constexpr Vector3 operator-(const Vector3& vector3) const noexcept
		{
			return Vector3(X - vector3.X, Y - vector3.Y, Z - vector3.Z);
		}

		constexpr Vector3 operator*(const Vector3& vector3) const noexcept
		{
			float x = X * vector3.X;
			float y = Y * vector3.Y;
			float z = Z * vector3.Z;

			return Vector3(x, y, z);
		}

		constexpr Vector3 operator/(const Vector3& vector3) const noexcept
		{
			float x = X / vector3.X;
			float y = Y / vector3.Y;
			float z = Z / vector3.Z;

			return Vector3(x, y, z);
		}

		constexpr Vector3 operator*(float f) const noexcept
		{
			return Vector3(X * f, Y * f, Z * f);
		}

		constexpr Vector3 operator/(float f) const noexcept
		{
			f = 1.0f / f;
			return Vector3(X * f, Y * f, Z * f);
		}

		Vector3& operator+=(const Vector3& vector3) noexcept
		{
			X += vector3.X;
			Y += vector3.Y;
			Z += vector3.Z;

			return *this;
		}

		Vector3& operator-=(const Vector3& vector3) noexcept
		{
			X -= vector3.X;
			Y -= vector3.Y;
			Z -= vector3.Z;

			return *this;
		}

		Vector3& operator+=(float f) noexcept
		{
			X += f;
			Y += f;
			Z += f;

			return *this;
		}

		Vector3& operator-=(float f) noexcept
		{
			X -= f;
			Y -= f;
			Z -= f;

			return *this;
		}

		Vector3& operator*=(float f) noexcept
		{
			X *= f;
			Y *= f;
			Z *= f;

			return *this;
		}

		Vector3& operator/=(float f) noexcept
		{
			f = 1.0f / f;

			X *= f;
			Y *= f;
			Z *= f;

			return *this;
		}

		constexpr bool operator==(const Vector3& vector3) const noexcept
		{
			return X == vector3.X && Y == vector3.Y && Z == vector3.Z;
		}

		constexpr bool operator!=(const Vector3& vector3) const noexcept
		{
			return !(vector3 == *this);
		}
	};

	// 4D Vector definition
	struct Vector4
	{
		float X{};
		float Y{};
		float Z{};
		float W{};

		constexpr Vector4() noexcept = default;
		constexpr Vector4(float x, float y, float z, float w) noexcept
			: X(x), Y(y), Z(z), W(w)
		{
		}

		constexpr explicit Vector4(const Vector3& vector3) noexcept
		{
			X = vector3.X;
			Y = vector3.Y;
			Z = vector3.Z;
			W = 0;
		}

		float& operator[](int i) noexcept
		{
			assert(i >= 0 && i < 4 && "Out of range: Vector4::Operator[]\n");
			return (&X)[i];
		}

		constexpr const float& operator[](int i) const noexcept
		{
			assert(i >= 0 && i < 4 && "Out of range: Vector4::Operator[]\n");
			return (&X)[i];
		}

		constexpr Vector4 operator-() const noexcept
		{
			return Vector4(-X, -Y, -Z, -W);
		}

		constexpr Vector4 operator+(const Vector4& vector4) const noexcept
		{
			return Vector4(X + vector4.X, Y + vector4.Y, Z + vector4.Z, W + vector4.W);
		}

		constexpr Vector4 operator-(const Vector4& vector4) const noexcept
		{
			return Vector4(X - vector4.X, Y - vector4.Y, Z - vector4.Z, W - vector4.W);
		}

		constexpr Vector4 operator*(const Vector4& vector4) const noexcept
		{
			float x = X * vector4.X;
			float y = Y * vector4.Y;
			float z = Z * vector4.Z;
			float w = W * vector4.W;

			return Vector4(x, y, z, w);
		}

		constexpr Vector4 operator/(const Vector4& vector4) const noexcept
		{
			float x = X / vector4.X;
			float y = Y / vector4.Y;
			float z = Z / vector4.Z;
			float w = W / vector4.W;

			return Vector4(x, y, z, w);
		}

		constexpr Vector4 operator*(float f) const noexcept
		{
			return Vector4(X * f, Y * f, Z * f, W * f);
		}

		constexpr Vector4 operator/(float f) const noexcept
		{
			f = 1.0f / f;
			return Vector4(X * f, Y * f, Z * f, W * f);
		}

		Vector4& operator+=(const Vector4& vector4) noexcept
		{
			X += vector4.X;
			Y += vector4.Y;
			Z += vector4.Z;
			W += vector4.W;

			return *this;
		}

		Vector4& operator-=(const Vector4& vector4) noexcept
		{
			X -= vector4.X;
			Y -= vector4.Y;
			Z -= vector4.Z;
			W -= vector4.W;

			return *this;
		}

		Vector4& operator+=(float f) noexcept
		{
			X += f;
			Y += f;
			Z += f;
			W += f;

			return *this;
		}

		Vector4& operator-=(float f) noexcept
		{
			X -= f;
			Y -= f;
			Z -= f;
			W -= f;

			return *this;
		}

		Vector4& operator*=(float f) noexcept
		{
			X *= f;
			Y *= f;
			Z *= f;
			W *= f;

			return *this;
		}

		Vector4& operator/=(float f) noexcept
		{
			f = 1.0f / f;
			X *= f;
			Y *= f;
			Z *= f;
			W *= f;

			return *this;
		}

		constexpr bool operator==(const Vector4& vector4) const noexcept
		{
			return X == vector4.X && Y == vector4.Y && Z == vector4.Z && W == vector4.W;
		}

		constexpr bool operator!=(const Vector4& vector4) const noexcept
		{
			return !(vector4 == *this);
		}
	};

	// 4x4 Matrix definition
	struct Matrix
	{
		constexpr Matrix() noexcept = default;
		constexpr explicit Matrix(float f) noexcept
		{
			mMatrix[0] = f;
			mMatrix[1] = 0;
			mMatrix[2] = 0;
			mMatrix[3] = 0;
			mMatrix[4] = 0;
			mMatrix[5] = f;
			mMatrix[6] = 0;
			mMatrix[7] = 0;
			mMatrix[8] = 0;
			mMatrix[9] = 0;
			mMatrix[10] = f;
			mMatrix[11] = 0;
			mMatrix[12] = 0;
			mMatrix[13] = 0;
			mMatrix[14] = 0;
			mMatrix[15] = f;
		}

		explicit Matrix(float* m) noexcept
		{
			memcpy(mMatrix, m, sizeof(float) * 16);
		}

		constexpr Matrix(float m00, float m01, float m02, float m03, float m04,
			float m05, float m06, float m07, float m08, float m09,
			float m10, float m11, float m12, float m13, float m14,
			float m15) noexcept
		{
			mMatrix[0] = m00;
			mMatrix[1] = m01;
			mMatrix[2] = m02;
			mMatrix[3] = m03;
			mMatrix[4] = m04;
			mMatrix[5] = m05;
			mMatrix[6] = m06;
			mMatrix[7] = m07;
			mMatrix[8] = m08;
			mMatrix[9] = m09;
			mMatrix[10] = m10;
			mMatrix[11] = m11;
			mMatrix[12] = m12;
			mMatrix[13] = m13;
			mMatrix[14] = m14;
			mMatrix[15] = m15;
		}

		float& operator[](unsigned int m) noexcept
		{
			assert(m >= 0 && m < 16 && "Out of range: Matrix::Operator[]\n");
			return mMatrix[m];
		}

		constexpr const float& operator[](unsigned int m) const noexcept
		{
			assert(m >= 0 && m < 16 && "Out of range: Matrix::Operator[]\n");
			return mMatrix[m];
		}

		float& operator()(unsigned int col, unsigned int row) noexcept
		{
			assert(col >= 0 && col < 4 && row >= 0 && row < 4 && "Out of range: Matrix::Operator[]\n");
			return mMatrix[4 * col + row];
		}

		constexpr const float& operator()(unsigned int col,
			unsigned int row) const noexcept
		{
			assert(col >= 0 && col < 4 && row >= 0 && row < 4 && "Out of range: Matrix::Operator[]\n");
			return mMatrix[4 * col + row];
		}

		constexpr Matrix operator*(const Matrix& matrix) const noexcept
		{
			float a00 = mMatrix[0];
			float a01 = mMatrix[1];
			float a02 = mMatrix[2];
			float a03 = mMatrix[3];
			float a04 = mMatrix[4];
			float a05 = mMatrix[5];
			float a06 = mMatrix[6];
			float a07 = mMatrix[7];
			float a08 = mMatrix[8];
			float a09 = mMatrix[9];
			float a10 = mMatrix[10];
			float a11 = mMatrix[11];
			float a12 = mMatrix[12];
			float a13 = mMatrix[13];
			float a14 = mMatrix[14];
			float a15 = mMatrix[15];

			float b00 = matrix[0];
			float b01 = matrix[1];
			float b02 = matrix[2];
			float b03 = matrix[3];
			float b04 = matrix[4];
			float b05 = matrix[5];
			float b06 = matrix[6];
			float b07 = matrix[7];
			float b08 = matrix[8];
			float b09 = matrix[9];
			float b10 = matrix[10];
			float b11 = matrix[11];
			float b12 = matrix[12];
			float b13 = matrix[13];
			float b14 = matrix[14];
			float b15 = matrix[15];

			float m00 = (a00 * b00) + (a04 * b01) + (a08 * b02) + (a12 * b03);
			float m01 = (a01 * b00) + (a05 * b01) + (a09 * b02) + (a13 * b03);
			float m02 = (a02 * b00) + (a06 * b01) + (a10 * b02) + (a14 * b03);
			float m03 = (a03 * b00) + (a07 * b01) + (a11 * b02) + (a15 * b03);
			float m04 = (a00 * b04) + (a04 * b05) + (a08 * b06) + (a12 * b07);
			float m05 = (a01 * b04) + (a05 * b05) + (a09 * b06) + (a13 * b07);
			float m06 = (a02 * b04) + (a06 * b05) + (a10 * b06) + (a14 * b07);
			float m07 = (a03 * b04) + (a07 * b05) + (a11 * b06) + (a15 * b07);
			float m08 = (a00 * b08) + (a04 * b09) + (a08 * b10) + (a12 * b11);
			float m09 = (a01 * b08) + (a05 * b09) + (a09 * b10) + (a13 * b11);
			float m10 = (a02 * b08) + (a06 * b09) + (a10 * b10) + (a14 * b11);
			float m11 = (a03 * b08) + (a07 * b09) + (a11 * b10) + (a15 * b11);
			float m12 = (a00 * b12) + (a04 * b13) + (a08 * b14) + (a12 * b15);
			float m13 = (a01 * b12) + (a05 * b13) + (a09 * b14) + (a13 * b15);
			float m14 = (a02 * b12) + (a06 * b13) + (a10 * b14) + (a14 * b15);
			float m15 = (a03 * b12) + (a07 * b13) + (a11 * b14) + (a15 * b15);

			return Matrix(m00, m01, m02, m03, m04, m05, m06, m07, m08, m09, m10, m11, m12,
				m13, m14, m15);
		}

		constexpr Matrix operator*(float f) const noexcept
		{
			float m00 = mMatrix[0] * f;
			float m01 = mMatrix[1] * f;
			float m02 = mMatrix[2] * f;
			float m03 = mMatrix[3] * f;
			float m04 = mMatrix[4] * f;
			float m05 = mMatrix[5] * f;
			float m06 = mMatrix[6] * f;
			float m07 = mMatrix[7] * f;
			float m08 = mMatrix[8] * f;
			float m09 = mMatrix[9] * f;
			float m10 = mMatrix[10] * f;
			float m11 = mMatrix[11] * f;
			float m12 = mMatrix[12] * f;
			float m13 = mMatrix[13] * f;
			float m14 = mMatrix[14] * f;
			float m15 = mMatrix[15] * f;

			return Matrix(m00, m01, m02, m03, m04, m05, m06, m07, m08, m09, m10, m11, m12,
				m13, m14, m15);
		}

		float* Data() noexcept
		{
			return mMatrix;
		}

		constexpr const float* Data() const noexcept
		{
			return mMatrix;
		}

	private:
		float mMatrix[16]{};
	};

	// Quaternion definition
	struct Quaternion
	{
		float X{};
		float Y{};
		float Z{};
		float W{};

		constexpr Quaternion() noexcept = default;
		constexpr Quaternion(float x, float y, float z, float w) noexcept
			: X(x), Y(y), Z(z), W(w)
		{
		}

		constexpr Quaternion(const Vector3& vector3, float w) noexcept
		{
			X = vector3.X;
			Y = vector3.Y;
			Z = vector3.Z;
			W = w;
		}

		float& operator[](int i) noexcept
		{
			assert(i >= 0 && i < 4 && "Out of range: Quaternion::Operator[]\n");
			return (&X)[i];
		}

		constexpr const float& operator[](int i) const noexcept
		{
			assert(i >= 0 && i < 4 && "Out of range: Quaternion::Operator[]\n");
			return (&X)[i];
		}

		constexpr Quaternion operator-() const noexcept
		{
			return Quaternion(-X, -Y, -Z, -W);
		}

		constexpr Quaternion operator+(const Quaternion& quaternion) const noexcept
		{
			return Quaternion(X + quaternion.X, Y + quaternion.Y, Z + quaternion.Z,
				W + quaternion.W);
		}

		constexpr Quaternion operator-(const Quaternion& quaternion) const noexcept
		{
			return Quaternion(X - quaternion.X, Y - quaternion.Y, Z - quaternion.Z,
				W - quaternion.W);
		}

		constexpr Quaternion operator*(const Quaternion& quaternion) const noexcept
		{
			float lhs_x = X;
			float lhs_y = Y;
			float lhs_z = Z;
			float lhs_w = W;

			float rhs_x = quaternion.X;
			float rhs_y = quaternion.Y;
			float rhs_z = quaternion.Z;
			float rhs_w = quaternion.W;

			float x =
				(lhs_x * rhs_w) + (lhs_w * rhs_x) + (lhs_y * rhs_z) - (lhs_z * rhs_y);
			float y =
				(lhs_y * rhs_w) + (lhs_w * rhs_y) + (lhs_z * rhs_x) - (lhs_x * rhs_z);
			float z =
				(lhs_z * rhs_w) + (lhs_w * rhs_z) + (lhs_x * rhs_y) - (lhs_y * rhs_x);
			float w =
				(lhs_w * rhs_w) - (lhs_x * rhs_x) - (lhs_y * rhs_y) - (lhs_z * rhs_z);

			return Quaternion(x, y, z, w);
		}

		constexpr Quaternion operator+(float f) const noexcept
		{
			return Quaternion(X + f, Y + f, Z + f, W + f);
		}

		constexpr Quaternion operator-(float f) const noexcept
		{
			return Quaternion(X - f, Y - f, Z - f, W - f);
		}

		Quaternion& operator+=(float f) noexcept
		{
			X += f;
			Y += f;
			Z += f;
			W += f;

			return *this;
		}

		Quaternion& operator-=(float f) noexcept
		{
			X -= f;
			Y -= f;
			Z -= f;
			W -= f;

			return *this;
		}

		Quaternion& operator*=(float f) noexcept
		{
			X *= f;
			Y *= f;
			Z *= f;
			W *= f;

			return *this;
		}

		Quaternion& operator/=(float f) noexcept
		{
			X /= f;
			Y /= f;
			Z /= f;
			W /= f;

			return *this;
		}

		constexpr Quaternion operator*(float f) const noexcept
		{
			return Quaternion(X * f, Y * f, Z * f, W * f);
		}

		constexpr Quaternion operator/(float f) const noexcept
		{
			f = 1.0f / f;
			return Quaternion(X * f, Y * f, Z * f, W * f);
		}

		constexpr bool operator==(const Quaternion& quaternion) const noexcept
		{
			return X == quaternion.X && Y == quaternion.Y && Z == quaternion.Z &&
				   W == quaternion.W;
		}

		constexpr bool operator!=(const Quaternion& quaternion) const noexcept
		{
			return !(quaternion == *this);
		}
	};

	// Plane definition
	struct Plane
	{
		Vector3 Normal{};
		float Offset{};

		Plane() = default;
		Plane(const Vector3& normal, float offset)
			: Normal(normal), Offset(offset) {}
	};

	// Ray definition
	struct Ray
	{
		Vector3 Origin{};
		Vector3 Direction{};

		constexpr Ray() noexcept = default;
		constexpr Ray(const Vector3& origin, const Vector3& direction) noexcept
			: Origin(origin), Direction(direction) {}

		constexpr Vector3 At(float t) noexcept
		{
			return Origin * (Direction * t);
		}
	};

	// 2D Vector
	constexpr Vector2 operator*(float f, const Vector2& vector2) noexcept
	{
		return Vector2(vector2.X * f, vector2.Y * f);
	}

	constexpr Vector2 operator/(float f, const Vector2& vector2) noexcept
	{
		f = 1.0f / f;
		return Vector2(vector2.X * f, vector2.Y * f);
	}

	inline float Distance(const Vector2& lhs, const Vector2& rhs) noexcept
	{
		float x = rhs.X - lhs.X;
		float y = rhs.Y - lhs.Y;

		return Sqrt((x * x) + (y * y));
	}

	constexpr float Dot(const Vector2& lhs, const Vector2& rhs) noexcept
	{
		return (lhs.X * rhs.X) + (lhs.Y * rhs.Y);
	}

	constexpr Vector2 Lerp(const Vector2& lhs, const Vector2& rhs, float t) noexcept
	{
		t = Clamp(0, 1.0f, t);

		float x = Lerp(lhs.X, rhs.X, t);
		float y = Lerp(lhs.Y, rhs.Y, t);

		return Vector2(x, y);
	}

	inline float Magnitude(const Vector2& vector2) noexcept
	{
		return std::sqrt((vector2.X * vector2.X) + (vector2.Y * vector2.Y));
	}

	constexpr float MagnitudeSquared(const Vector2& vector2) noexcept
	{
		return (vector2.X * vector2.X) + (vector2.Y * vector2.Y);
	}

	inline Vector2 Normalized(const Vector2& vector2) noexcept
	{
		return vector2 / Magnitude(vector2);
	}

	constexpr Vector2 Reflect(const Vector2& vector2, const Vector2& normal) noexcept
	{
		float dot_product = Dot(vector2, normal);
		float x = vector2.X - ((2.0f * normal.X) * dot_product);
		float y = vector2.Y - ((2.0f * normal.Y) * dot_product);

		return Vector2(x, y);
	}

	inline Vector2 Rotate(const Vector2& vector2, float r) noexcept
	{
		float cos = Cos(r);
		float sin = Sin(r);
		float x = (vector2.X * cos) - (vector2.Y * sin);
		float y = (vector2.X * sin) + (vector2.Y * cos);

		return Vector2(x, y);
	}

	constexpr Vector2 Scale(const Vector2& vector2, float f) noexcept
	{
		return vector2 * f;
	}

	constexpr Vector2 Transform(const Vector2& vector2, const Matrix& matrix) noexcept
	{
		float x = (matrix[0] * vector2.X) + (matrix[4] * vector2.Y) + (matrix[8]) +
				  (matrix[12]);
		float y = (matrix[1] * vector2.X) + (matrix[5] * vector2.Y) + (matrix[9]) +
				  (matrix[13]);

		return Vector2(x, y);
	}

	// 3D Vector
	constexpr Vector3 operator*(float f, const Vector3& vector3) noexcept
	{
		return Vector3(vector3.X * f, vector3.Y * f, vector3.Z * f);
	}

	constexpr Vector3 operator/(float f, const Vector3& vector3) noexcept
	{
		f = 1.0f / f;
		return Vector3(vector3.X * f, vector3.Y * f, vector3.Z * f);
	}

	constexpr Vector3 Cross(const Vector3& lhs, const Vector3& rhs) noexcept
	{
		float x = (lhs.Y * rhs.Z) - (lhs.Z * rhs.Y);
		float y = (lhs.Z * rhs.X) - (lhs.X * rhs.Z);
		float z = (lhs.X * rhs.Y) - (lhs.Y * rhs.X);

		return Vector3(x, y, z);
	}

	inline float Distance(const Vector3& lhs, const Vector3& rhs) noexcept
	{
		float x = rhs.X - lhs.X;
		float y = rhs.Y - lhs.Y;
		float z = rhs.Z - lhs.Z;

		return Sqrt((x * x) + (y * y) + (z * z));
	}

	constexpr float Dot(const Vector3& lhs, const Vector3& rhs) noexcept
	{
		return (lhs.X * rhs.X) + (lhs.Y * rhs.Y) + (lhs.Z * rhs.Z);
	}

	constexpr Vector3 Lerp(const Vector3& lhs, const Vector3& rhs, float t) noexcept
	{
		t = Clamp(0, 1.0f, t);

		float x = Lerp(lhs.X, rhs.X, t);
		float y = Lerp(lhs.Y, rhs.Y, t);
		float z = Lerp(lhs.Z, rhs.Z, t);

		return Vector3(x, y, z);
	}

	float Magnitude(const Vector3& vector3) noexcept
	{
		return Sqrt((vector3.X * vector3.X) + (vector3.Y * vector3.Y) +
					(vector3.Z * vector3.Z));
	}

	constexpr float MagnitudeSquared(const Vector3& vector3) noexcept
	{
		return (vector3.X * vector3.X) + (vector3.Y * vector3.Y) +
			   (vector3.Z * vector3.Z);
	}

	Vector3 Normalized(const Vector3& vector3) noexcept
	{
		return vector3 / Magnitude(vector3);
	}

	constexpr Vector3 Reflect(const Vector3& vector3, const Vector3& normal) noexcept
	{
		float dot_product = Dot(vector3, normal);
		float x = vector3.X - ((2.0f * normal.X) * dot_product);
		float y = vector3.Y - ((2.0f * normal.Y) * dot_product);
		float z = vector3.Z - ((2.0f * normal.Z) * dot_product);

		return Vector3(x, y, z);
	}

	inline Vector3 Rotate(const Vector3& vector3, const Quaternion& quaternion) noexcept
	{
		Vector3 qv = Vector3(quaternion.X, quaternion.Y, quaternion.Z);

		return vector3 * ((quaternion.W * quaternion.W) - Dot(qv, qv)) +
			   (2.0f * qv * Dot(qv, vector3)) +
			   (2.0f * quaternion.W * Cross(qv, vector3));
	}

	constexpr Vector3 Scale(const Vector3& vector3, float f) noexcept
	{
		return vector3 * f;
	}

	inline Vector3 Transform(const Vector3& vector3, const Matrix& matrix) noexcept
	{
		return Vector3((matrix[0] * vector3.X) + (matrix[4] * vector3.Y) +
						   (matrix[8] * vector3.Z) + matrix[12],
			(matrix[1] * vector3.X) + (matrix[5] * vector3.Y) +
				(matrix[9] * vector3.Z) + matrix[13],
			(matrix[2] * vector3.X) + (matrix[6] * vector3.Y) +
				(matrix[10] * vector3.Z) + matrix[14]);
	}

	// 4D Vector
	inline float Distance(const Vector4& lhs, const Vector4& rhs) noexcept
	{
		float x = rhs.X - lhs.X;
		float y = rhs.Y - lhs.Y;
		float z = rhs.Z - lhs.Z;
		float w = rhs.W - lhs.W;

		return Sqrt((x * x) + (y * y) + (z * z) + (w * w));
	}

	constexpr float Dot(const Vector4& lhs, const Vector4& rhs) noexcept
	{
		return (lhs.X * rhs.X) + (lhs.Y * rhs.Y) + (lhs.Z * rhs.Z) + (lhs.W * rhs.W);
	}

	constexpr Vector4 Lerp(const Vector4& lhs, const Vector4& rhs, float t) noexcept
	{
		t = Clamp(0, 1.0f, t);

		float x = Lerp(lhs.X, rhs.X, t);
		float y = Lerp(lhs.Y, rhs.Y, t);
		float z = Lerp(lhs.Z, rhs.Z, t);
		float w = Lerp(lhs.W, rhs.W, t);

		return Vector4(x, y, z, w);
	}

	inline float Magnitude(const Vector4& vector4) noexcept
	{
		return Sqrt((vector4.X * vector4.X) + (vector4.Y * vector4.Y) +
					(vector4.Z * vector4.Z) + (vector4.W * vector4.W));
	}

	constexpr float MagnitudeSquared(const Vector4& vector4) noexcept
	{
		return (vector4.X * vector4.X) + (vector4.Y * vector4.Y) +
			   (vector4.Z * vector4.Z) + (vector4.W * vector4.W);
	}

	inline Vector4 Normalized(const Vector4& vector4) noexcept
	{
		return vector4 / Magnitude(vector4);
	}

	constexpr Vector4 Reflect(const Vector4& vector4, const Vector4& normal) noexcept
	{
		float dot_product = Dot(vector4, normal);
		float x = vector4.X - ((2.0f * normal.X) * dot_product);
		float y = vector4.Y - ((2.0f * normal.Y) * dot_product);
		float z = vector4.Z - ((2.0f * normal.Z) * dot_product);
		float w = vector4.W - ((2.0f * normal.W) * dot_product);

		return Vector4(x, y, z, w);
	}

	constexpr Vector4 operator*(float f, const Vector4& vector4) noexcept
	{
		return Vector4(vector4.X * f, vector4.Y * f, vector4.Z * f, vector4.W * f);
	}

	constexpr Vector4 operator/(float f, const Vector4& vector4) noexcept
	{
		f = 1.0f / f;
		return Vector4(vector4.X * f, vector4.Y * f, vector4.Z * f, vector4.W * f);
	}

	// 4x4 Matrix
	static const Matrix IdentityMatrix = Matrix(1.0f, 0, 0, 0, 0, 1.0f, 0, 0, 0, 0, 1.0f, 0, 0, 0, 0, 1.0f);

	constexpr float Determinant(const Matrix& matrix) noexcept
	{
		float m00 = matrix[0];
		float m01 = matrix[1];
		float m02 = matrix[2];
		float m03 = matrix[3];
		float m04 = matrix[4];
		float m05 = matrix[5];
		float m06 = matrix[6];
		float m07 = matrix[7];
		float m08 = matrix[8];
		float m09 = matrix[9];
		float m10 = matrix[10];
		float m11 = matrix[11];
		float m12 = matrix[12];
		float m13 = matrix[13];
		float m14 = matrix[14];
		float m15 = matrix[15];

		float a = (m12 * m09 * m06 * m03) - (m08 * m13 * m06 * m03) -
				  (m12 * m05 * m10 * m03) + (m04 * m13 * m10 * m03);
		float b = (m08 * m05 * m14 * m03) - (m04 * m09 * m14 * m03) -
				  (m12 * m09 * m02 * m07) + (m08 * m13 * m02 * m07);
		float c = (m12 * m01 * m10 * m07) - (m00 * m13 * m10 * m07) -
				  (m08 * m01 * m14 * m07) + (m00 * m09 * m14 * m07);
		float d = (m12 * m05 * m02 * m11) - (m04 * m13 * m02 * m11) -
				  (m12 * m01 * m06 * m11) + (m00 * m13 * m06 * m11);
		float e = (m04 * m01 * m14 * m11) - (m00 * m05 * m14 * m11) -
				  (m08 * m05 * m02 * m15) + (m04 * m09 * m02 * m15);
		float f = (m08 * m01 * m06 * m15) - (m00 * m09 * m06 * m15) -
				  (m04 * m01 * m10 * m15) + (m00 * m05 * m10 * m15);

		return a + b + c + d + e + f;
	}

	constexpr Matrix Inverse(const Matrix& matrix) noexcept
	{
		float d = Determinant(matrix);

		if (d == 0)
		{
			return Matrix();
		}

		float id = 1.0f / d;

		float m00 = matrix[0];
		float m01 = matrix[1];
		float m02 = matrix[2];
		float m03 = matrix[3];
		float m04 = matrix[4];
		float m05 = matrix[5];
		float m06 = matrix[6];
		float m07 = matrix[7];
		float m08 = matrix[8];
		float m09 = matrix[9];
		float m10 = matrix[10];
		float m11 = matrix[11];
		float m12 = matrix[12];
		float m13 = matrix[13];
		float m14 = matrix[14];
		float m15 = matrix[15];

		float im00 = ((m05 * m10 * m15) + (m09 * m14 * m07) + (m13 * m06 * m11) -
						 (m05 * m14 * m11) - (m09 * m06 * m15) - (m13 * m10 * m07)) *
					 id;

		float im01 = ((m13 * m10 * m03) - (m09 * m14 * m03) - (m13 * m02 * m11) +
						 (m01 * m14 * m11) + (m09 * m02 * m15) - (m01 * m10 * m15)) *
					 id;

		float im02 = ((m05 * m14 * m03) - (m13 * m06 * m03) + (m13 * m02 * m07) -
						 (m01 * m14 * m07) - (m05 * m02 * m15) + (m01 * m06 * m15)) *
					 id;

		float im03 = ((m09 * m06 * m03) - (m05 * m10 * m03) - (m09 * m02 * m07) +
						 (m01 * m10 * m07) + (m05 * m02 * m11) - (m01 * m06 * m11)) *
					 id;

		float im04 = ((m12 * m10 * m07) - (m08 * m14 * m07) - (m12 * m06 * m11) +
						 (m04 * m14 * m11) + (m08 * m06 * m15) - (m04 * m10 * m15)) *
					 id;

		float im05 = ((m08 * m14 * m03) - (m12 * m10 * m03) + (m12 * m02 * m11) -
						 (m00 * m14 * m11) - (m08 * m02 * m15) + (m00 * m10 * m15)) *
					 id;

		float im06 = ((m12 * m06 * m03) - (m04 * m14 * m03) - (m12 * m02 * m07) +
						 (m00 * m14 * m07) + (m04 * m02 * m15) - (m00 * m06 * m15)) *
					 id;

		float im07 = ((m04 * m10 * m03) - (m08 * m06 * m03) + (m08 * m02 * m07) -
						 (m00 * m10 * m07) - (m04 * m02 * m11) + (m00 * m06 * m11)) *
					 id;

		float im08 = ((m08 * m13 * m07) - (m12 * m09 * m07) + (m12 * m05 * m11) -
						 (m04 * m13 * m11) - (m08 * m05 * m15) + (m04 * m09 * m15)) *
					 id;

		float im09 = ((m12 * m09 * m03) - (m08 * m13 * m03) - (m12 * m01 * m11) +
						 (m00 * m13 * m11) + (m08 * m01 * m15) - (m00 * m09 * m15)) *
					 id;

		float im10 = ((m04 * m13 * m03) - (m12 * m05 * m03) + (m12 * m01 * m07) -
						 (m00 * m13 * m07) - (m04 * m01 * m15) + (m00 * m05 * m15)) *
					 id;

		float im11 = ((m08 * m05 * m03) - (m04 * m09 * m03) - (m08 * m01 * m07) +
						 (m00 * m09 * m07) + (m04 * m01 * m11) - (m00 * m05 * m11)) *
					 id;

		float im12 = ((m12 * m09 * m06) - (m08 * m13 * m06) - (m12 * m05 * m10) +
						 (m04 * m13 * m10) + (m08 * m05 * m14) - (m04 * m09 * m14)) *
					 id;

		float im13 = ((m08 * m13 * m02) - (m12 * m09 * m02) + (m12 * m01 * m10) -
						 (m00 * m13 * m10) - (m08 * m01 * m14) + (m00 * m09 * m14)) *
					 id;

		float im14 = ((m12 * m05 * m02) - (m04 * m13 * m02) - (m12 * m01 * m06) +
						 (m00 * m13 * m06) + (m04 * m01 * m14) - (m00 * m05 * m14)) *
					 id;

		float im15 = ((m04 * m09 * m02) - (m08 * m05 * m02) + (m08 * m01 * m06) -
						 (m00 * m09 * m06) - (m04 * m01 * m10) + (m00 * m05 * m10)) *
					 id;

		return Matrix(im00, im01, im02, im03, im04, im05, im06, im07, im08, im09,
			im10, im11, im12, im13, im14, im15);
	}

	constexpr Matrix Lerp(const Matrix& lhs, const Matrix& rhs, float t) noexcept
	{
		float m00 = Lerp(lhs[0], rhs[0], t);
		float m01 = Lerp(lhs[1], rhs[1], t);
		float m02 = Lerp(lhs[2], rhs[2], t);
		float m03 = Lerp(lhs[3], rhs[3], t);
		float m04 = Lerp(lhs[4], rhs[4], t);
		float m05 = Lerp(lhs[5], rhs[5], t);
		float m06 = Lerp(lhs[6], rhs[6], t);
		float m07 = Lerp(lhs[7], rhs[7], t);
		float m08 = Lerp(lhs[8], rhs[8], t);
		float m09 = Lerp(lhs[9], rhs[9], t);
		float m10 = Lerp(lhs[10], rhs[10], t);
		float m11 = Lerp(lhs[11], rhs[11], t);
		float m12 = Lerp(lhs[12], rhs[12], t);
		float m13 = Lerp(lhs[13], rhs[13], t);
		float m14 = Lerp(lhs[14], rhs[14], t);
		float m15 = Lerp(lhs[15], rhs[15], t);

		return Matrix(m00, m01, m02, m03, m04, m05, m06, m07, m08, m09, m10, m11, m12,
			m13, m14, m15);
	}

	inline Matrix LookAt(const Vector3& position, const Vector3& target,
		const Vector3& up) noexcept
	{
		Vector3 z = sg::Normalized(position - target);
		Vector3 x = sg::Normalized(sg::Cross(up, z));
		Vector3 y = sg::Cross(z, x);

		Matrix matrix = IdentityMatrix;

		matrix[0] = x.X;
		matrix[1] = y.X;
		matrix[2] = z.X;
		matrix[3] = 0;
		matrix[4] = x.Y;
		matrix[5] = y.Y;
		matrix[6] = z.Y;
		matrix[7] = 0;
		matrix[8] = x.Z;
		matrix[9] = y.Z;
		matrix[10] = z.Z;
		matrix[11] = 0;
		matrix[12] = -Dot(x, position);
		matrix[13] = -Dot(y, position);
		matrix[14] = -Dot(z, position);
		matrix[15] = 1.0f;

		return matrix;
	}

	constexpr Matrix MakeFrustum(float l, float r, float b, float t, float n,
		float f) noexcept
	{
		float rl = 1.0f / (r - l);
		float tb = 1.0f / (t - b);
		float fn = 1.0f / (f - n);

		float m00 = (n * 2.0f) * rl;
		float m05 = (n * 2.0f) * tb;
		float m08 = (r + l) * rl;
		float m09 = (t + b) * tb;
		float m10 = -(f + n) * fn;
		float m14 = -(f * n * 2.0f) * fn;

		return Matrix(m00, 0, 0, 0, 0, m05, 0, 0, m08, m09, m10, -1.0f, 0, 0, m14, 0);
	}

	constexpr Matrix MakeOrthographic(float l, float r, float t, float b, float n,
		float f) noexcept
	{
		float w = 1.0f / (r - l);
		float h = 1.0f / (b - t);
		float p = 1.0f / (f - n);

		float x = (r + l) * w;
		float y = (t + b) * h;
		float z = (f + n) * p;

		return Matrix(2.0f * w, 0, 0, -x, 0, 2.0f * h, 0, -y, 0, 0, p, -z, 0, 0, 0,
			1.0f);
	}

	inline Matrix MakePerspective(float fovy, float a, float n, float f) noexcept
	{
		float t = n * Tan(fovy * 0.5f);
		float r = t * a;
		return MakeFrustum(-r, r, -t, t, n, f);
	}

	inline Matrix MakeRotationFromAxisAngle(float x, float y, float z, float r) noexcept
	{
		float magnitude = Sqrt((x * x) + (y * y) + (z * z));

		if ((magnitude != 1.0f) && (magnitude != 0.0f))
		{
			magnitude = 1.0f / magnitude;
			x *= magnitude;
			y *= magnitude;
			z *= magnitude;
		}

		float cos = Cos(r);
		float sin = Sin(r);

		float t = 1.0f - cos;

		float tx = t * x;
		float ty = t * y;
		float tz = t * z;

		float m00 = (tx * x) + cos;
		float m01 = (tx * y) + (sin * z);
		float m02 = (tx * z) - (sin * y);
		float m04 = (tx * y) - (sin * z);
		float m05 = (ty * y) + cos;
		float m06 = (ty * z) + (sin * x);
		float m08 = (tx * z) + (sin * y);
		float m09 = (ty * z) - (sin * x);
		float m10 = (tz * z) + cos;

		return Matrix(m00, m01, m02, 0, m04, m05, m06, 0, m08, m09, m10, 0, 0, 0, 0,
			1.0f);
	}

	inline Matrix MakeRotationFromAxisAngle(const Vector3& axis, float r) noexcept
	{
		return MakeRotationFromAxisAngle(axis.X, axis.Y, axis.Z, r);
	}

	inline Matrix MakeRotationFromAxisAngle(const Matrix& matrix,
		const Vector3& axis, float r) noexcept
	{
		return matrix * MakeRotationFromAxisAngle(axis, r);
	}

	inline Matrix MakeRotationFromEuler(float roll, float pitch, float yaw) noexcept
	{
		float a = Cos(roll);
		float b = Cos(pitch);
		float c = Cos(yaw);
		float d = Sin(roll);
		float e = Sin(pitch);
		float f = Sin(yaw);

		float ac = a * c;
		float af = a * f;
		float dc = d * c;
		float df = d * f;

		float m00 = b * c;
		float m01 = -b * f;
		float m02 = e;
		float m04 = af + dc * e;
		float m05 = ac - df * e;
		float m06 = -d * b;
		float m08 = df - ac * e;
		float m09 = dc + af * e;
		float m10 = a * b;

		return Matrix(m00, m01, m02, 0, m04, m05, m06, 0, m08, m09, m10, 0, 0, 0, 0,
			1.0f);
	}

	inline Matrix MakeRotationFromEuler(const Vector3& euler) noexcept
	{
		return MakeRotationFromEuler(euler.X, euler.Y, euler.Z);
	}

	inline Matrix MakeRotationX(float r) noexcept
	{
		Matrix matrix4 = IdentityMatrix;
		float cos = Cos(r);
		float sin = Sin(r);

		matrix4[5] = cos;
		matrix4[6] = -sin;
		matrix4[9] = sin;
		matrix4[10] = cos;

		return matrix4;
	}

	inline Matrix MakeRotationY(float r) noexcept
	{
		Matrix matrix4 = IdentityMatrix;
		float cos = Cos(r);
		float sin = Sin(r);

		matrix4[0] = cos;
		matrix4[2] = sin;
		matrix4[8] = -sin;
		matrix4[10] = cos;

		return matrix4;
	}

	inline Matrix MakeRotationZ(float r) noexcept
	{
		Matrix matrix4 = IdentityMatrix;
		float cos = Cos(r);
		float sin = Sin(r);

		matrix4[0] = cos;
		matrix4[1] = -sin;
		matrix4[4] = sin;
		matrix4[5] = cos;

		return matrix4;
	}

	constexpr Matrix MakeScale(float x, float y, float z) noexcept
	{
		return Matrix(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1.0f);
	}

	constexpr Matrix MakeScale(const Vector3& vector3) noexcept
	{
		return Matrix(vector3.X, 0, 0, 0, 0, vector3.Y, 0, 0, 0, 0, vector3.Z, 0, 0,
			0, 0, 1.0f);
	}

	constexpr Matrix MakeShear(float x, float y, float z) noexcept
	{
		return Matrix(
			1.0f, x, x, 0,
			y, 1.0f, y, 0,
			z, z, 1.0f, 0,
			0, 0, 0, 1.0f);
	}

	constexpr Matrix MakeShear(const Vector3& vector3) noexcept
	{
		return MakeShear(vector3.X, vector3.Y, vector3.Z);
	}

	constexpr Matrix MakeTranslationMatrix(float x, float y, float z) noexcept
	{
		return Matrix(1.0f, 0, 0, 0, 0, 1.0f, 0, 0, 0, 0, 1.0f, 0, x, y, z, 1.0f);
	}

	constexpr Matrix MakeTranslationMatrix(const Vector3& vector3) noexcept
	{
		return MakeTranslationMatrix(vector3.X, vector3.Y, vector3.Z);
	}

	constexpr Matrix Transposed(const Matrix& matrix) noexcept
	{
		return Matrix(matrix[0], matrix[1], matrix[2], matrix[3], matrix[4],
			matrix[5], matrix[6], matrix[7], matrix[8], matrix[9],
			matrix[10], matrix[11], matrix[12], matrix[13], matrix[14],
			matrix[15]);
	}

	constexpr Matrix operator*(float f, const Matrix& matrix4) noexcept
	{
		return matrix4 * f;
	}

	// Quaternion
	constexpr Quaternion operator*(float f, const Quaternion& quaternion) noexcept
	{
		return Quaternion(quaternion.X * f, quaternion.Y * f, quaternion.Z * f,
			quaternion.W * f);
	}

	constexpr Quaternion operator/(float f, const Quaternion& quaternion) noexcept
	{
		f = 1.0f / f;
		return Quaternion(quaternion.X * f, quaternion.Y * f, quaternion.Z * f,
			quaternion.W * f);
	}

	static const Quaternion IdentityQuaternion = Quaternion(0, 0, 0, 1.0f);

	constexpr float Dot(const Quaternion& lhs, const Quaternion& rhs) noexcept
	{
		return (lhs.X * rhs.X) + (lhs.Y * rhs.Y) + (lhs.Z * rhs.Z) + (lhs.W * rhs.W);
	}

	inline float AngleBetween(const Quaternion& lhs, const Quaternion& rhs) noexcept
	{
		float clamped = Abs(Clamp(-1, 1, Dot(lhs, rhs)));
		return 2 * Acos(clamped);
	}

	inline float Magnitude(const Quaternion& quaternion) noexcept
	{
		return Sqrt((quaternion.X * quaternion.X) + (quaternion.Y * quaternion.Y) +
					(quaternion.Z * quaternion.Z) + (quaternion.W * quaternion.W));
	}

	constexpr float MagnitudeSquared(const Quaternion& quaternion) noexcept
	{
		return (quaternion.X * quaternion.X) + (quaternion.Y * quaternion.Y) +
			   (quaternion.Z * quaternion.Z) + (quaternion.W * quaternion.W);
	}

	inline Quaternion Normalized(const Quaternion& quaternion) noexcept
	{
		float f = 1.0f / Magnitude(quaternion);
		return Quaternion(quaternion.X * f, quaternion.X * f, quaternion.X * f,
			quaternion.X * f);
	}

	constexpr Quaternion Lerp(const Quaternion& lhs, const Quaternion& rhs, float t) noexcept
	{
		t = Clamp(0, 1.0f, t);

		float x = Lerp(lhs.X, rhs.X, t);
		float y = Lerp(lhs.Y, rhs.Y, t);
		float z = Lerp(lhs.Z, rhs.Z, t);
		float w = Lerp(lhs.W, rhs.W, t);

		return Quaternion(x, y, z, w);
	}

	inline Quaternion NLerp(const Quaternion& lhs, const Quaternion& rhs, float t) noexcept
	{
		t = Clamp(0, 1.0f, t);
		return Normalized(Lerp(lhs, rhs, t));
	}

	inline Quaternion Slerp(const Quaternion& lhs, const Quaternion& rhs, float t) noexcept
	{
		t = Clamp(0, 1.0f, t);
		float dot = Dot(lhs, rhs);

		if (Abs(dot) >= 1.0f)
		{
			return lhs;
		}

		if (dot > 0.95f)
		{
			return NLerp(lhs, rhs, t);
		}

		float half_angle = Acos(dot);
		float sin_half_angle = 1.0f / Sin(half_angle);

		if (Abs(sin_half_angle) < Epsilon)
		{
			return (lhs * 0.5) + (rhs * 0.5);
		}

		float sin0 = Sin((1.0f - t) * half_angle) / sin_half_angle;
		float sin1 = Sin(t * half_angle) / sin_half_angle;

		return (lhs * sin0) + (rhs * sin1);
	}

	inline Quaternion MakeQuaternionFromAxisAngle(float x, float y, float z,
		float angle) noexcept
	{
		float half_angle = angle * 0.5f;
		float sin = Sin(half_angle);
		float cos = Cos(half_angle);

		float qx = x * sin;
		float qy = y * sin;
		float qz = z * sin;
		float qw = cos;

		return Quaternion(qx, qy, qz, qw);
	}

	inline Quaternion MakeQuaternionFromAxisAngle(const Vector3& axis,
		float angle) noexcept
	{
		return MakeQuaternionFromAxisAngle(axis.X, axis.Y, axis.Z, angle);
	}

	inline Quaternion MakeQuaternionFromEuler(float yaw, float pitch, float roll) noexcept
	{
		yaw *= 0.5;
		pitch *= 0.5;
		roll *= 0.5;

		float cos_yaw = Cos(yaw);
		float sin_yaw = Sin(yaw);
		float cos_pitch = Cos(pitch);
		float sin_pitch = Sin(pitch);
		float cos_roll = Cos(roll);
		float sin_roll = Sin(roll);

		float qx =
			(cos_yaw * sin_pitch * cos_roll) + (sin_yaw * cos_pitch * sin_roll);
		float qy =
			(sin_yaw * cos_pitch * cos_roll) - (cos_yaw * sin_pitch * sin_roll);
		float qz =
			(cos_yaw * cos_pitch * sin_roll) - (sin_yaw * sin_pitch * cos_roll);
		float qw =
			(cos_yaw * cos_pitch * cos_roll) + (sin_yaw * sin_pitch * sin_roll);

		return Quaternion(qx, qy, qz, qw);
	}

	inline Quaternion MakeQuaternionFromEuler(const Vector3& euler) noexcept
	{
		return MakeQuaternionFromEuler(euler.X, euler.Y, euler.Z);
	}

	inline Quaternion MakeQuaternionFromMatrix(const Matrix& matrix) noexcept
	{
		float m00 = matrix[0];
		float m10 = matrix[1];
		float m20 = matrix[2];
		float m01 = matrix[4];
		float m11 = matrix[5];
		float m21 = matrix[6];
		float m02 = matrix[8];
		float m12 = matrix[9];
		float m22 = matrix[10];

		float trace = m00 + m11 + m22;
		float s{};
		float qx{};
		float qy{};
		float qz{};
		float qw{};

		if (trace > 0)
		{
			qw = Sqrt(trace + 1.0f) * 0.5f;

			s = 0.25f / qw;

			qx = (m21 - m12) * s;
			qy = (m02 - m20) * s;
			qz = (m10 - m01) * s;
		}
		else if (m00 > m11 && m00 > m22)
		{
			qx = Sqrt(1.0f + m00 - m11 - m22);

			s = 0.25f / qx;

			qw = (m21 - m12) * s;
			qy = (m01 + m10) * s;
			qz = (m02 + m20) * s;
		}
		else if (m11 > m22)
		{
			qy = Sqrt(1.0f + m11 - m00 - m22);

			s = 0.25f / qy;

			qw = (m02 - m20) * s;
			qx = (m01 + m10) * s;
			qz = (m12 + m21) * s;
		}
		else
		{
			qz = Sqrt(1.0f + m22 - m00 - m11);

			s = 0.25f / qz;

			qx = (m02 + m20) * s;
			qy = (m12 + m21) * s;
			qw = (m10 - m01) * s;
		}
		return Quaternion(qx, qy, qz, qw);
	}

	constexpr Matrix ToMatrix(const Quaternion& quaternion) noexcept
	{
		float xx = quaternion.X * quaternion.X;
		float yy = quaternion.Y * quaternion.Y;
		float zz = quaternion.Z * quaternion.Z;
		float wx = quaternion.W * quaternion.X;
		float wy = quaternion.W * quaternion.Y;
		float wz = quaternion.W * quaternion.Z;
		float xy = quaternion.X * quaternion.Y;
		float xz = quaternion.X * quaternion.Z;
		float yz = quaternion.Y * quaternion.Z;

		float m00 = 1.0F - (2.0F * (yy + zz));
		float m10 = 2.0F * (xy - wz);
		float m20 = 2.0F * (xz + wy);
		float m01 = 2.0F * (xy + wz);
		float m11 = 1.0F - (2.0F * (xx + zz));
		float m21 = 2.0F * (yz - wx);
		float m02 = 2.0F * (xz - wy);
		float m12 = 2.0F * (yz + wx);
		float m22 = 1.0F - (2.0F * (xx + yy));

		return Matrix(m00, m10, m20, 0, m01, m11, m21, 0, m02, m12, m22, 0, 0, 0, 0,
			1.0f);
	}

	// Ray
	inline float Intersects(const Ray& ray, const Plane& plane) noexcept
	{
		float d = Dot(plane.Normal, ray.Direction);

		if (Abs(d) <= Epsilon)
		{
			float n = Dot(plane.Normal, ray.Origin) + plane.Offset;
			return -(n / d);
		}

		return 0;
	}

} // namespace sg

#endif // SGMATH_H
