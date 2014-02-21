#ifndef __BGK3D15V_HPP_INCLUDED__
#define __BGK3D15V_HPP_INCLUDED__

#include<cstdio>
#include"field.hpp"

/**
 * BGK3D15V
 */
template<std::size_t N,std::size_t M,std::size_t L>
class BGK3D15V{

	public:
		static const std::size_t A = 15;	//速度ベクトルの個数
		static const std::size_t DIM = 3;	//次元
		static const int C[A][DIM];			//速度ベクトル
		static const double E[A];			//feqで使うやつ
		static const std::size_t BACK[A];	//反対方向のindex

	private:
		Field<double,N,M,L> u,v,w;			//速度
		Field<double,N,M,L> rho;			//密度
		Field<double,N,M,L> f[A+A];				//分布関数
		const double tau_inv;				//時定数の逆数
		std::size_t step;					//計算ステップ数

	public:
		BGK3D15V(double t)					//コンストラクタ
			: tau_inv(1./t), step(0){
				initial();
			}
		virtual ~BGK3D15V(){}

		void output(FILE* fp);					//ファイル出力
		virtual void evolution(std::size_t n);	//計算をnステップ進める

	protected:
		virtual double feq(std::size_t a,
				std::size_t x,std::size_t y,std::size_t z) const;
		//局所平衡分布
		virtual void collision();					//衝突
		void streaming(Field<double,N,M,L> ff[]);	//移動
		virtual void macros();						//マクロな値を計算する
		virtual void initial();						//初期条件
		virtual void boundary();					//境界条件
};



/******************** 静的メンバ変数の初期化 ********************/

template<std::size_t N,std::size_t M,std::size_t L>
const int BGK3D15V<N,M,L>::C[A][DIM] = { 
		{ 0, 0, 0},		//0
		//最近接格子
		{ 1, 0, 0},		//1
		{ 0, 1, 0},		//2
		{-1, 0, 0},		//3
		{ 0,-1, 0},		//4
		{ 0, 0, 1},		//5
		{ 0, 0,-1},		//6
		//次近接格子
		{ 1, 1, 1},		//7
		{-1, 1, 1},		//8
		{-1,-1, 1},		//9
		{ 1,-1, 1},		//10
		{ 1, 1,-1},		//11
		{-1, 1,-1},		//12
		{-1,-1,-1},		//13
		{ 1,-1,-1},		//14
	};

	template<std::size_t N,std::size_t M,std::size_t L>
		const double BGK3D15V<N,M,L>::E[A] = { 2./9.,
			1./9., 1./9., 1./9., 1./9., 1./9., 1./9.,
			1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};

	template<std::size_t N,std::size_t M,std::size_t L>
		const std::size_t BGK3D15V<N,M,L>::BACK[A] = { 0, 
			3, 4, 1, 2, 6, 5, 
			13, 14, 11, 12, 9, 10, 7, 8
		}; 




	/******************** メソッドの実装 ********************/

	/**
	 * ファイルへの出力
	 *
	 * x y z u v w ρ の順のcsvを出力
	 * x,y はindexなので整数
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::output(FILE* fp){
			for(std::size_t x=0; x<N; x++){
				for(std::size_t y=0; y<M; y++){
					for(std::size_t z=0; z<L; z++){
						std::fprintf(fp, "%lu\t%lu\t%lu\t%e\t%e\t%e\t%e\n", 
								x, y, z, 
								u(x,y,z), v(x,y,z), w(x,y,z), rho(x,y,z));
					}
				}
				std::fprintf(fp,"\n\n");
			}
		}


	/**
	 * 時間発展
	 *
	 * main関数でコイツを呼ぶ
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::evolution(std::size_t n){
			while(n--){
				//衝突
				collision();

				//移動
				streaming(f);

				//境界条件
				boundary();

				//マクロな値の計算
				macros();

				step++;
			}
		}


	/**
	 * 局所平衡分布 (BGK3D15Vモデル)
	 *
	 * @param a 速度のindex
	 * @prama x xのindex
	 * @param y yのindex
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline double BGK3D15V<N,M,L>::feq(std::size_t a,
				std::size_t x, std::size_t y, std::size_t z) const{

			const double u_ = u(x,y,z);
			const double v_ = v(x,y,z);
			const double w_ = w(x,y,z);
			const double cdotu 
				= C[a][0]*u(x,y,z) + C[a][1]*v(x,y,z) + C[a][2]*w(x,y,z);
			return E[a]*rho(x,y,z)*(
					1. + 
					+ 3.* cdotu
					+ 9./2.* cdotu * cdotu
					- 3./2.*(u_*u_ + v_*v_ + w_*w_)
					);
		}


	/**
	 * 粒子の衝突
	 *
	 * 格子ボルツマン方程式を計算する
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::collision(){
			std::size_t x;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(x=0; x<N; x++){
				for(std::size_t y=0; y<M; y++){
					for(std::size_t z=0; z<L; z++){
						for(std::size_t a=0; a<A; a++){
							f[a+A](x,y,z) 
								= (1.-tau_inv)*f[a](x,y,z)
								+ tau_inv*feq(a,x,y,z);
						}
					}
				}
			}
		}


	/**
	 * 粒子の移動
	 *
	 * 壁際を特別扱いしてることに注意
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::streaming(Field<double,N,M,L> ff[A]){
			std::size_t x;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(x=0; x<N; x++){
				for(std::size_t y=0; y<M; y++){
					for(std::size_t z=0; z<L; z++){
						for(std::size_t a=0; a<A; a++){
							ff[a](x+C[a][0],y+C[a][1],z+C[a][2]) = ff[a+A](x,y,z);
						}
					}
				}
			}
		}


	/**
	 * 巨視的変数の計算
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::macros(){
			std::size_t x;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(x=0;x<N;x++){
				for(std::size_t y=0;y<M; y++){
					for(std::size_t z=0; z<L; z++){
						double rho_ = 0;	//密度の一時的な変数
						double ru_ = 0;		//運動量の一時的な変数
						double rv_ = 0;		//同上
						double rw_ = 0;		//同上

						for(std::size_t a=0;a<A;a++){
							//密度の計算
							rho_ += f[a](x,y,z);

							//運動量の計算
							ru_ += C[a][0]*f[a](x,y,z);
							rv_ += C[a][1]*f[a](x,y,z);
							rw_ += C[a][2]*f[a](x,y,z);
						}
						//代入
						rho(x,y,z) = rho_;
						u(x,y,z) = ru_/rho_;
						v(x,y,z) = rv_/rho_;
						w(x,y,z) = rw_/rho_;
					}
				}
			}
		}


	/**
	 * 境界条件
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::boundary(){
			std::size_t x,y,z;

			//周期境界 with 圧力勾配
			const static double dp = 0.01;	//圧力差
			const static std::size_t x0 = 0;
			const static std::size_t xL = N-1;
			for(y=0; y<M; y++){
				for(z=0; z<L; z++){
					const double c = dp - 1./3.*(
							  f[0+A](x0,y,z) - f[0+A](xL,y,z)
							+ f[2+A](x0,y,z) - f[2+A](xL,y,z)
							+ f[4+A](x0,y,z) - f[4+A](xL,y,z)
							+ f[5+A](x0,y,z) - f[5+A](xL,y,z)
							+ f[6+A](x0,y,z) - f[6+A](xL,y,z)
							);
					//左端
					f[ 1](x0,y,z) = f[ 1](xL,y,z) + c;
					f[ 7](x0,y,z) = f[ 7](xL,y,z) + c/8.;
					f[10](x0,y,z) = f[10](xL,y,z) + c/8.;
					f[11](x0,y,z) = f[11](xL,y,z) + c/8.;
					f[14](x0,y,z) = f[14](xL,y,z) + c/8.;

					//右端
					f[ 3](xL,y,z) = f[ 3](x0,y,z) - c;
					f[ 8](xL,y,z) = f[ 8](x0,y,z) - c/8.;
					f[ 9](xL,y,z) = f[ 9](x0,y,z) - c/8.;
					f[12](xL,y,z) = f[12](x0,y,z) - c/8.;
					f[13](xL,y,z) = f[13](x0,y,z) - c/8.;
				}
			}

			//固定壁 at 北
			y=M-1;
			for(x=0; x<N; x++){
				for(z=0; z<L; z++){
					f[4](x,y,z) = f[BACK[4]+A](x,y,z);
					f[9](x,y,z) = f[BACK[9]+A](x,y,z);
					f[10](x,y,z) = f[BACK[10]+A](x,y,z);
					f[13](x,y,z) = f[BACK[13]+A](x,y,z);
					f[14](x,y,z) = f[BACK[14]+A](x,y,z);
				}
			}

			//固定壁 at 南
			y=0;
			for(x=0; x<N; x++){
				for(z=0; z<L; z++){
					f[2](x,y,z) = f[BACK[2]+A](x,y,z);
					f[7](x,y,z) = f[BACK[7]+A](x,y,z);
					f[8](x,y,z) = f[BACK[8]+A](x,y,z);
					f[11](x,y,z) = f[BACK[11]+A](x,y,z);
					f[12](x,y,z) = f[BACK[12]+A](x,y,z);
				}
			}

			//固定壁 at 上
			z = L-1;
			for(x=0; x<N; x++){
				for(y=0; y<M; y++){
					f[6](x,y,z) = f[BACK[6]+A](x,y,z);
					f[11](x,y,z) = f[BACK[11]+A](x,y,z);
					f[12](x,y,z) = f[BACK[12]+A](x,y,z);
					f[13](x,y,z) = f[BACK[13]+A](x,y,z);
					f[14](x,y,z) = f[BACK[14]+A](x,y,z);
				}
			}

			//固定壁 at 下
			z = 0;
			for(x=0; x<N; x++){
				for(y=0; y<M; y++){
					f[5](x,y,z) = f[BACK[5]+A](x,y,z);
					f[7](x,y,z) = f[BACK[7]+A](x,y,z);
					f[8](x,y,z) = f[BACK[8]+A](x,y,z);
					f[9](x,y,z) = f[BACK[9]+A](x,y,z);
					f[10](x,y,z) = f[BACK[10]+A](x,y,z);
				}
			}
		}


	/**
	 * 初期条件
	 */
	template<std::size_t N,std::size_t M,std::size_t L>
		inline void BGK3D15V<N,M,L>::initial(){
			//各点で密度が１になるようにする
			for(std::size_t x=0; x<N; x++){
				for(std::size_t y=0; y<M; y++){
					for(std::size_t z=0; z<M; z++){
						for(std::size_t a=0; a<A; a++){
							f[a](x,y,z) = 1./A;
						}
					}
				}
			}
			//マクロな値を計算する
			macros();
		}
#endif //__BGK3D15V_HPP_INCLUDED__
