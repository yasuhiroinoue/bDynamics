/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef QUATERNION_H
#define QUATERNION_H

#include <stdio.h>
#include <cmath>
#include <iostream>
#include "vec.h"

using namespace std;

template <class TEMPLATE>
class _quat{
	public:
	TEMPLATE n;
	TEMPLATE x;
	TEMPLATE y;
	TEMPLATE z;
	_quat(void);
	_quat(const TEMPLATE&, const TEMPLATE&, const TEMPLATE&, const TEMPLATE& );
	_quat(const TEMPLATE&, const TEMPLATE&, const TEMPLATE&);
	_quat(const TEMPLATE&, const _vec<TEMPLATE>&);
	_quat(const _quat<TEMPLATE>&);
	
	
	
	//オペレーター
	void IN(const TEMPLATE& , const TEMPLATE& , const TEMPLATE&, const TEMPLATE& );
	void IN(const TEMPLATE& , const TEMPLATE& , const TEMPLATE& );
	void IN(const TEMPLATE&, const _vec<TEMPLATE>&);
	void IN(const _quat<TEMPLATE>&);
	_quat<int> icast(void) const;
	_quat<TEMPLATE>& operator = (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator += (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator -= (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator *= (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator *= (const TEMPLATE&);
	_quat<TEMPLATE>& operator /= (const TEMPLATE&);
	
	
	TEMPLATE norm(void) const;
	TEMPLATE angl(void) const;
	_vec<TEMPLATE> vec(void) const;
	_vec<TEMPLATE> axis(void) const;
	_quat<TEMPLATE> qrot(const _quat<TEMPLATE>&, const _quat<TEMPLATE>& );
	_vec<TEMPLATE>  vrot(const _quat<TEMPLATE>&, const  _vec<TEMPLATE>& );
	
	
	_quat<TEMPLATE> qmake( const TEMPLATE&, const TEMPLATE&, const TEMPLATE& );
};






template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(void):n(0),x(0),y(0),z(0){}

template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(const TEMPLATE& d, const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c):n(d),x(a),y(b),z(c){}

template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	TEMPLATE ca = std::cos(a/2.0);		TEMPLATE sa = std::sin(a/2.0);
	TEMPLATE cb = std::cos(b/2.0);		TEMPLATE sb = std::sin(b/2.0);
	TEMPLATE cc = std::cos(c/2.0);		TEMPLATE sc = std::sin(c/2.0);
	
	n = (ca*cb*cc + sa*sb*sc);
	x = (sa*cb*cc - ca*sb*sc);
	y = (ca*sb*cc + sa*cb*sc);
	z = (ca*cb*sc - sa*sb*cc);
}

template <class T> inline _quat<T>::_quat(const T& angle, const _vec<T>& v){
	n = cos(angle*0.5);
	T tmp = sin(angle*0.5);
	x = tmp * v.x;
	y = tmp * v.y;
	z = tmp * v.z;
}

template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(const _quat<TEMPLATE>& q){
	n = q.n;
	x = q.x;
	y = q.y;
	z = q.z;
}

//オペレータ
template <class TEMPLATE> inline void _quat<TEMPLATE>::IN(const TEMPLATE& d, const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	
	n = d; x = a; y = b; z = c;
}
template <class TEMPLATE> inline void _quat<TEMPLATE>::IN(const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	TEMPLATE ca = std::cos(a/2.0);		TEMPLATE sa = std::sin(a/2.0);
	TEMPLATE cb = std::cos(b/2.0);		TEMPLATE sb = std::sin(b/2.0);
	TEMPLATE cc = std::cos(c/2.0);		TEMPLATE sc = std::sin(c/2.0);
	
	n = (ca*cb*cc + sa*sb*sc);
	x = (sa*cb*cc - ca*sb*sc);
	y = (ca*sb*cc + sa*cb*sc);
	z = (ca*cb*sc - sa*sb*cc);
}
template <class T> inline void _quat<T>::IN(const T& angle, const _vec<T>& v){
	n = std::cos(angle*0.5);
	T tmp = std::sin(angle*0.5);
	x = tmp * v.x;
	y = tmp * v.y;
	z = tmp * v.z;
}

//template <class TEMPLATE> inline void _quat<TEMPLATE>::IN(const _quat<TEMPLATE>& v){
//	n = v.s; x = v.x; y = v.y; z = v.z;
//}
template <class TEMPLATE> inline void _quat<TEMPLATE>::IN(const _quat<TEMPLATE>& v){
	n = v.n; x = v.x; y = v.y; z = v.z;
}//deji make a modification//2010/05/25

//クォータニオン対スカラー

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator *= ( const TEMPLATE& a){
	n *= a; x *= a; y *= a; z *= a;
	return *this;
}

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator /= ( const TEMPLATE& a){
	n /= a; x /= a; y /= a; z /= a;
	return *this;
}


//クォータニオン対クォータニオン
template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator =  ( const _quat<TEMPLATE>& v ){
  n = v.n;
  x = v.x;
  y = v.y;
  z = v.z;
  return *this;
}

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator +=  ( const _quat<TEMPLATE>& v ){
  n += v.n;
  x += v.x;
  y += v.y;
  z += v.z;
  return *this;
}

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator -=  ( const _quat<TEMPLATE>& v ){
  n -= v.n;
  x -= v.x;
  y -= v.y;
  z -= v.z;
  return *this;
}


//乗算
template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator *= ( const _quat<TEMPLATE>& q2){
	_quat<TEMPLATE> q1(n,x,y,z);
	
	n = q1.n * q2.n - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	x = q1.n * q2.x + q2.n * q1.x + q1.y * q2.z - q1.z * q2.y;
	y = q1.n * q2.y + q2.n * q1.y + q1.z * q2.x - q1.x * q2.z;
	z = q1.n * q2.z + q2.n * q1.z + q1.x * q2.y - q1.y * q2.x;
	
	return *this;
}


//スカラ対クォータニオン
template <class TEMPLATE> inline _quat<TEMPLATE> operator*( const TEMPLATE& a, const _quat<TEMPLATE>& q  ){
	_quat<TEMPLATE> quat(a*q.n, a*q.x, a*q.y, a*q.z);
	return quat;
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator*( const _quat<TEMPLATE>& q, const TEMPLATE& a  ){
	_quat<TEMPLATE> quat(a*q.n, q.x*a, q.y*a, q.z*a);
	return quat;
}

template <class TEMPLATE> inline _quat<TEMPLATE> operator/(const _quat<TEMPLATE>& q, const TEMPLATE& a  ){
	_quat<TEMPLATE> quat(q.n/a, q.x/a, q.y/a, q.z/a);
	return quat;
}


//クォータニオン対ベクトル
template <class TEMPLATE> inline _quat<TEMPLATE> operator*(const _quat<TEMPLATE>& q1, const _vec<TEMPLATE>& q2){
	_quat<TEMPLATE> quat;
	
	quat.n = - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	quat.x = q1.n * q2.x  + q1.y * q2.z - q1.z * q2.y;
	quat.y = q1.n * q2.y  + q1.z * q2.x - q1.x * q2.z;
	quat.z = q1.n * q2.z  + q1.x * q2.y - q1.y * q2.x;
	
	return quat;
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator*(const _vec<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	_quat<TEMPLATE> quat;
	quat.n =  - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	quat.x = q2.n * q1.x + q1.y * q2.z - q1.z * q2.y;
	quat.y = q2.n * q1.y + q1.z * q2.x - q1.x * q2.z;
	quat.z = q2.n * q1.z + q1.x * q2.y - q1.y * q2.x;
	
	return quat;
}

//クォータニオン対クォータニオン
template <class TEMPLATE> inline _quat<TEMPLATE> operator*(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	_quat<TEMPLATE> quat;
	quat.n = q1.n * q2.n - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	quat.x = q1.n * q2.x + q2.n * q1.x + q1.y * q2.z - q1.z * q2.y;
	quat.y = q1.n * q2.y + q2.n * q1.y + q1.z * q2.x - q1.x * q2.z;
	quat.z = q1.n * q2.z + q2.n * q1.z + q1.x * q2.y - q1.y * q2.x;
	
	return quat;
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator+(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	
	return _quat<TEMPLATE>(q1.n + q2.n, q1.x+q2.x, q1.y+q2.y, q1.z + q2.z);
}

template <class TEMPLATE> inline _quat<TEMPLATE> operator-(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	return _quat<TEMPLATE>(q1.n - q2.n, q1.x - q2.x, q1.y - q2.y, q1.z - q2.z);
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator~(const _quat<TEMPLATE>& q1){
	return _quat<TEMPLATE>(q1.n, -q1.x, -q1.y, -q1.z);
}


template <class TEMPLATE> inline TEMPLATE _quat<TEMPLATE>::norm(void) const{
	TEMPLATE s = n*n + x*x + y*y + z*z;
	return sqrt(s);
}


template <class TEMPLATE> inline TEMPLATE _quat<TEMPLATE>::angl(void) const{
	return (2.0*acos(n));
}

template <class TEMPLATE> inline _vec<TEMPLATE> _quat<TEMPLATE>::vec(void) const{
	return _vec<TEMPLATE> (x,y,z);
}

template <class TEMPLATE> inline _vec<TEMPLATE> _quat<TEMPLATE>::axis(void) const{
	_vec<TEMPLATE> v = this->vec();
	TEMPLATE s = v.norm();
	
	if( s > 0 ) return v/s;
}


template <class TEMPLATE> inline _quat<TEMPLATE> _quat<TEMPLATE>::qrot(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	return q1*q2*(~q1);
}

template <class TEMPLATE> inline _vec<TEMPLATE> _quat<TEMPLATE>::vrot(const _quat<TEMPLATE>& q1, const _vec<TEMPLATE>& v2){
	_quat<TEMPLATE> t = q1*v2*(~q1);
	return t.vec();
	
}

template <class TEMPLATE> inline _quat<TEMPLATE> _quat<TEMPLATE>::qmake( const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	TEMPLATE ca = cos(a/2.0);
	TEMPLATE sa = sin(a/2.0);
	
	TEMPLATE cb = cos(b/2.0);
	TEMPLATE sb = sin(b/2.0);
	
	TEMPLATE cc = cos(c/2.0);
	TEMPLATE sc = sin(c/2.0);
	
	_quat<TEMPLATE> q((ca*cb*cc + sa*sb*sc),  (sa*cb*cc - ca*sb*sc), (ca*sb*cc + sa*cb*sc), (ca*cb*sc - sa*sb*cc) );
	return q;
}

namespace _quaternion{
	template <class T> inline _quat<T> rot(const _quat<T>& q1, const _quat<T>& q2){
		return (q1*q2*(~q1));
	}
	template <class T> inline _vec<T> rot(const _quat<T>& q1, const _vec<T>& v2){
		return ((q1*v2*(~q1)).vec());
	}
	template <class T> inline _quat<T> make( const T& a, const T& b, const T& c){//θψφ
		T ca = std::cos(a/2.0);		T sa = std::sin(a/2.0);
		T cb = std::cos(b/2.0);		T sb = std::sin(b/2.0);
		T cc = std::cos(c/2.0);		T sc = std::sin(c/2.0);
		
		return _quat<T>((ca*cb*cc + sa*sb*sc),  (sa*cb*cc - ca*sb*sc), (ca*sb*cc + sa*cb*sc), (ca*cb*sc - sa*sb*cc) );
	}
}

template <class T> inline std::ostream& operator<< (std::ostream& stream,const _quat<T> quat){
	return stream<<quat.n<<'\t'<<quat.x<<'\t'<<quat.y<<'\t'<<quat.z;
}



#endif

/*#ifndef QUATERNION_H
#define QUATERNION_H

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "vec.h"

using namespace std;

template <class TEMPLATE>
class _quat{
	public:
	TEMPLATE n;
	TEMPLATE x;
	TEMPLATE y;
	TEMPLATE z;
	_quat(void);
	_quat(const TEMPLATE&, const TEMPLATE&, const TEMPLATE&, const TEMPLATE& );
	_quat(const _quat<TEMPLATE>&);
	
	
	
	//オペレーター
	void IN(const TEMPLATE& , const TEMPLATE& , const TEMPLATE&, const TEMPLATE& );
	void IN(const _quat<TEMPLATE>&);
	_quat<int> icast(void) const;
	_quat<TEMPLATE>& operator = (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator += (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator -= (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator *= (const _quat<TEMPLATE>&);
	_quat<TEMPLATE>& operator *= (const TEMPLATE&);
	_quat<TEMPLATE>& operator /= (const TEMPLATE&);
	
	
	TEMPLATE norm(void) const;
	TEMPLATE angl(void) const;
	_vec<TEMPLATE> vec(void) const;
	_vec<TEMPLATE> axis(void) const;
	_quat<TEMPLATE> qrot(const _quat<TEMPLATE>, const _quat<TEMPLATE> );
	_vec<TEMPLATE>  vrot(const _quat<TEMPLATE>, const  _vec<TEMPLATE> );
	
	
	_quat<TEMPLATE> qmake( const TEMPLATE, const TEMPLATE, const TEMPLATE );
};






template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(void){
	n = 0;
	x = 0;
	y = 0;
	z = 0;
}

template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(const TEMPLATE& d, const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	n = d;
	x = a;
	y = b;
	z = c;
}

template <class TEMPLATE> inline _quat<TEMPLATE>::_quat(const _quat<TEMPLATE>& q){
	n = q.n;
	x = q.x;
	y = q.y;
	z = q.z;
}

//オペレータ
template <class TEMPLATE> inline void _quat<TEMPLATE>::IN(const TEMPLATE& d, const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	
	n = d; x = a; y = b; z = c;
}

template <class TEMPLATE> inline void _quat<TEMPLATE>::IN(const _quat<TEMPLATE>& v){
	n = v.s; x = v.x; y = v.y; z = v.z;
}

//クォータニオン対スカラー

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator *= ( const TEMPLATE& a){
	n *= a; x *= a; y *= a; z *= a;
	return *this;
}

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator /= ( const TEMPLATE& a){
	n /= a; x /= a; y /= a; z /= a;
	return *this;
}


//クォータニオン対クォータニオン
template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator =  ( const _quat<TEMPLATE>& v ){
  n = v.n;
  x = v.x;
  y = v.y;
  z = v.z;
  return *this;
}

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator +=  ( const _quat<TEMPLATE>& v ){
  n += v.n;
  x += v.x;
  y += v.y;
  z += v.z;
  return *this;
}

template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator -=  ( const _quat<TEMPLATE>& v ){
  n -= v.n;
  x -= v.x;
  y -= v.y;
  z -= v.z;
  return *this;
}


//乗算
template <class TEMPLATE> inline _quat<TEMPLATE>& _quat<TEMPLATE>::operator *= ( const _quat<TEMPLATE>& q2){
	_quat<TEMPLATE> q1(n,x,y,z);
	
	n = q1.n * q2.n - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	x = q1.n * q2.x + q2.n * q1.x + q1.y * q2.z - q1.z * q2.y;
	y = q1.n * q2.y + q2.n * q1.y + q1.z * q2.x - q1.x * q2.z;
	z = q1.n * q2.z + q2.n * q1.z + q1.x * q2.y - q1.y * q2.x;
	
	return *this;
}


//スカラ対クォータニオン
template <class TEMPLATE> inline _quat<TEMPLATE> operator*( const TEMPLATE a, const _quat<TEMPLATE> q  ){
	_quat<TEMPLATE> quat(a*q.n, a*q.x, a*q.y, a*q.z);
	return quat;
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator*( const _quat<TEMPLATE> q, const TEMPLATE a  ){
	_quat<TEMPLATE> quat(a*q.n, q.x*a, q.y*a, q.z*a);
	return quat;
}

template <class TEMPLATE> inline _quat<TEMPLATE> operator/(const _quat<TEMPLATE> q, const TEMPLATE a  ){
	_quat<TEMPLATE> quat(q.n/a, q.x/a, q.y/a, q.z/a);
	return quat;
}


//クォータニオン対ベクトル
template <class TEMPLATE> inline _quat<TEMPLATE> operator*(const _quat<TEMPLATE>& q1, const _vec<TEMPLATE>& q2){
	_quat<TEMPLATE> quat;
	
	quat.n = - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	quat.x = q1.n * q2.x  + q1.y * q2.z - q1.z * q2.y;
	quat.y = q1.n * q2.y  + q1.z * q2.x - q1.x * q2.z;
	quat.z = q1.n * q2.z  + q1.x * q2.y - q1.y * q2.x;
	
	return quat;
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator*(const _vec<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	_quat<TEMPLATE> quat;
	quat.n =  - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	quat.x = q2.n * q1.x + q1.y * q2.z - q1.z * q2.y;
	quat.y = q2.n * q1.y + q1.z * q2.x - q1.x * q2.z;
	quat.z = q2.n * q1.z + q1.x * q2.y - q1.y * q2.x;
	
	return quat;
}

//クォータニオン対クォータニオン
template <class TEMPLATE> inline _quat<TEMPLATE> operator*(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	_quat<TEMPLATE> quat;
	quat.n = q1.n * q2.n - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
	quat.x = q1.n * q2.x + q2.n * q1.x + q1.y * q2.z - q1.z * q2.y;
	quat.y = q1.n * q2.y + q2.n * q1.y + q1.z * q2.x - q1.x * q2.z;
	quat.z = q1.n * q2.z + q2.n * q1.z + q1.x * q2.y - q1.y * q2.x;
	
	return quat;
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator+(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	
	return _quat<TEMPLATE>(q1.n + q2.n, q1.x+q2.x, q1.y+q2.y, q1.z + q2.z);
}

template <class TEMPLATE> inline _quat<TEMPLATE> operator-(const _quat<TEMPLATE>& q1, const _quat<TEMPLATE>& q2){
	return _quat<TEMPLATE>(q1.n - q2.n, q1.x - q2.x, q1.y - q2.y, q1.z - q2.z);
}


template <class TEMPLATE> inline _quat<TEMPLATE> operator~(const _quat<TEMPLATE> q1){
	return _quat<TEMPLATE>(q1.n, -q1.x, -q1.y, -q1.z);
}


template <class TEMPLATE> inline TEMPLATE _quat<TEMPLATE>::norm(void) const{
	TEMPLATE s = n*n + x*x + y*y + z*z;
	return sqrt(s);
}


template <class TEMPLATE> inline TEMPLATE _quat<TEMPLATE>::angl(void) const{
	return (2.0*acos(n));
}

template <class TEMPLATE> inline _vec<TEMPLATE> _quat<TEMPLATE>::vec(void) const{
	return _vec<TEMPLATE> (x,y,z);
}

template <class TEMPLATE> inline _vec<TEMPLATE> _quat<TEMPLATE>::axis(void) const{
	_vec<TEMPLATE> v = this->vector();
	TEMPLATE s = v.norm();
	
	if( s > 0 ) return v/s;
}


template <class TEMPLATE> inline _quat<TEMPLATE> _quat<TEMPLATE>::qrot(const _quat<TEMPLATE> q1, const _quat<TEMPLATE> q2){
	return q1*q2*(~q1);
}

template <class TEMPLATE> inline _vec<TEMPLATE> _quat<TEMPLATE>::vrot(const _quat<TEMPLATE> q1, const _vec<TEMPLATE> v2){
	_quat<TEMPLATE> t = q1*v2*(~q1);
	return t.vec();
	
}

template <class TEMPLATE> inline _quat<TEMPLATE> _quat<TEMPLATE>::qmake( const TEMPLATE a, const TEMPLATE b, const TEMPLATE c){
	TEMPLATE ca = cos(a/2.0);
	TEMPLATE sa = sin(a/2.0);
	
	TEMPLATE cb = cos(b/2.0);
	TEMPLATE sb = sin(b/2.0);
	
	TEMPLATE cc = cos(c/2.0);
	TEMPLATE sc = sin(c/2.0);
	
	_quat<TEMPLATE> q((ca*cb*cc + sa*sb*sc),  (sa*cb*cc - ca*sb*sc), (ca*sb*cc + sa*cb*sc), (ca*cb*sc - sa*sb*cc) );
	return q;
}


#endif

*/
