#ifndef __FIELD_HPP_INCLUDED__
#define __FIELD_HPP_INCLUDED__


#include<valarray>
#include<cmath>
#include<cstring>


/**
 * スカラー場
 *
 * 壁にめり込んでもいいようにちょっと余計にメモリを確保している
 *
 * T : 型
 * N,M,L : x,y,z方向の格子数
 */
template<typename T, std::size_t N,std::size_t M, std::size_t L=1>
class Field : public std::valarray<T> {
	private:
		static constexpr size_t W = 1;	//壁の分の大きさ

	public:
		static constexpr size_t NX = N;
		static constexpr size_t NY = M;
		static constexpr size_t NZ = L;

		/**
		 * コンストラクタ
		 */
		Field(){
			this->resize((NX+W+W)*(NY+W+W)*(NZ+W+W));
		}
		~Field(){}

		/**
		 * 整数引数
		 * 値の参照を返す
		 */
		T& operator() (std::size_t x,std::size_t y,std::size_t z=0){
			return (*this)[W+x+(y+W)*(NX+W+W)+(z+W)*(NX+W+W)*(NY+W+W)];
		}
		const T& operator() (std::size_t x,std::size_t y,std::size_t z=0) const{
			return (*this)[W+x+(y+W)*(NX+W+W)+(z+W)*(NX+W+W)*(NY+W+W)];
		}
};


#endif //__FIELD_HPP_INCLUDED__
