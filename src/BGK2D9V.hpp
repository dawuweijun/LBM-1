#ifndef __BGK2D9V_H_INCLUDED__
#define __BGK2D9V_H_INCLUDED__

#include<cstdio>
#include"field.hpp"



/**
 * BGK2D9Vモデル
 */
template<std::size_t N,std::size_t M>
class BGK2D9V{

	public:
		static const std::size_t A = 9;		//速度ベクトルの個数
		static const std::size_t DIM = 2;	//次元
		static const int C[A][DIM];			//速度ベクトル
		static const double E[A];			//feqで使うやつ
		static const std::size_t BACK[A];	//反対方向のindex

	private:
		Field<double,N,M,1> u,v;			//速度
		Field<double,N,M,1> rho;			//密度
		Field<double,N,M,1> f[A+A];			//分布関数
		const double tau_inv;				//時定数の逆数
		std::size_t step;					//計算ステップ数

	public:
		BGK2D9V(double t)					//コンストラクタ
			: tau_inv(1./t), step(0){
			initial();
		}
		virtual ~BGK2D9V(){}
		
		void output(FILE* fp);					//ファイル出力
		virtual void evolution(std::size_t n);	//計算をnステップ進める
		
	protected:
		virtual double feq(std::size_t a,
				std::size_t x,std::size_t y) const;	//局所平衡分布
		virtual void collision();					//衝突
		void streaming(Field<double,N,M,1> ff[]);	//移動
		virtual void macros();						//マクロな値を計算する
		virtual void initial();						//初期条件
		virtual void boundary();					//境界条件
};



/******************** 静的メンバ変数の初期化 ********************/

template<std::size_t N,std::size_t M>
const int BGK2D9V<N,M>::C[A][DIM] = {
	{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1},
};

template<std::size_t N,std::size_t M>
const double BGK2D9V<N,M>::E[A] = {
	4./9., 1./9., 1./9., 1./9., 1./9.,1./36.,1./36.,1./36.,1./36.
};

template<std::size_t N,std::size_t M>
const std::size_t BGK2D9V<N,M>::BACK[A] = {
	0, 3, 4, 1, 2, 7, 8, 5, 6
};




/******************** メソッドの実装 ********************/

/**
 * ファイルへの出力
 *
 * x y u v ρ の順のcsvを出力
 * x,y はindexなので整数
 */
template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::output(FILE* fp){
	for(std::size_t x=0; x<N; x++){
		for(std::size_t y=0; y<M; y++){
			std::fprintf(fp,"%d\t%d\t%e\t%e\t%e\n", 
					x, y, u(x,y), v(x,y), rho(x,y));
		}
	}
}


/**
 * 時間発展
 *
 * main関数でコイツを呼ぶ
 */
template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::evolution(std::size_t n){
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
 * 局所平衡分布 (BGK2D9Vモデル)
 *
 * @param a 速度のindex
 * @prama x xのindex
 * @param y yのindex
 */
template<std::size_t N,std::size_t M>
inline double BGK2D9V<N,M>::feq(std::size_t a,std::size_t x,std::size_t y) const{
	const double cdotu = C[a][0]*u(x,y) + C[a][1]*v(x,y);
	return E[a]*rho(x,y)*(
			1. + 
			+ 3.*(C[a][0]*u(x,y) + C[a][1]*v(x,y))
			+ 9./2.* cdotu * cdotu
			- 3./2.*(u(x,y)*u(x,y)+v(x,y)*v(x,y))
			);
}


/**
 * 粒子の衝突
 *
 * 格子ボルツマン方程式を計算する
 */
template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::collision(){
	//double sum=0;
	for(std::size_t x=0; x<N; x++){
		for(std::size_t y=0; y<M; y++){
			for(std::size_t a=0; a<A; a++){
				f[a+A](x,y) = (1.-tau_inv)*f[a](x,y) + tau_inv*feq(a,x,y);
	//			sum += -tau_inv*(f[a](x,y) - feq(a,x,y));
			}
		}
	}
}


/**
 * 粒子の移動
 *
 * 壁際を特別扱いしてることに注意
 */
template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::streaming(Field<double,N,M,1> ff[A]){
	std::size_t x,y,a;

	//バルクからの移動
	for(x=0; x<N; x++){
		for(y=0; y<M; y++){
			for(a=0; a<A; a++){
				ff[a](x+C[a][0],y+C[a][1]) = ff[a+A](x,y);
			}
		}
	}
}

template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::macros(){
	for(std::size_t x=0;x<N;x++){
		for(std::size_t y=0;y<M;y++){
			double rho_ = 0;	//密度のダミー変数
			double ru_ = 0;		//運動量のダミー変数
			double rv_ = 0;		//運動量のダミー変数

			for(std::size_t a=0;a<A;a++){
				//密度の計算
				rho_ += f[a](x,y);

				//運動量の計算
				ru_ += C[a][0]*f[a](x,y);
				rv_ += C[a][1]*f[a](x,y);
			}
			//代入
			rho(x,y) = rho_;
			u(x,y) = ru_/rho_;
			v(x,y) = rv_/rho_;
		}
	}
}


/**
 * 境界条件
 */
template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::boundary(){
	std::size_t x,y,a;

	//周期境界 with 圧力勾配
	const static double dp = 0.01;	//圧力差
	const static std::size_t x0 = 0;
	const static std::size_t xL = N-1;
	for(y=0; y<M; y++){
		const double c = dp - 1./3.*(
	 		  f[0+A](x0,y) - f[0+A](xL,y)
			+ f[2+A](x0,y) - f[2+A](xL,y)
			+ f[4+A](x0,y) - f[4+A](xL,y)
			);
		//左端
		f[1](x0,y) = f[1+A](xL,y) + c;
		f[5](x0,y) = f[5+A](xL,y) + c/4.;
		f[8](x0,y) = f[8+A](xL,y) + c/4.;

		//右端
		f[3](xL,y) = f[3+A](x0,y) - c;
		f[6](xL,y) = f[6+A](x0,y) - c/4.;
		f[7](xL,y) = f[7+A](x0,y) - c/4.;
	}

	//固定壁 at 上
	y=M-1;
	for(x=0; x<N; x++){
		f[4](x,y) = f[BACK[4]+A](x,y);
		f[7](x,y) = f[BACK[7]+A](x,y);
		f[8](x,y) = f[BACK[8]+A](x,y);
	}

	//固定壁 at 下
	y=0;
	for(x=0; x<N; x++){
		f[2](x,y) = f[BACK[2]+A](x,y);
		f[5](x,y) = f[BACK[5]+A](x,y);
		f[6](x,y) = f[BACK[6]+A](x,y);
	}
}


/**
 * 初期条件
 */
template<std::size_t N,std::size_t M>
inline void BGK2D9V<N,M>::initial(){
	//各点で密度が１になるようにする
	for(std::size_t x=0; x<N; x++){
		for(std::size_t y=0; y<M; y++){
			for(std::size_t a=0; a<A; a++){
				f[a](x,y) = 1./A;
			}
		}
	}
	//マクロな値を計算する
	macros();
}


#endif //__BGK2D9V_H_INCLUDED__
