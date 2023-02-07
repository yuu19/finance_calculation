#include <iostream>
#include <string>
#include <cmath>
using namespace std;
//二項ツリー
// r:株価 sigma:ボラティリティ T:満期 K:行使価格 S0:初期時点での原資産の価格 N:分割数
double r, sigma, T, K, S0, N;


// op = 1 コール op = 0　プット
double binary_tree(double r, double sigma, double T, double K, double S0, double N, int op = 1)

{
  //株価行列
  vector<vector<double>> S(N+1, vector<double>(N+1, 0));
  //オプション価格行列
  vector<vector<double>> S(N+1, vector<double>(N+1, 0));

  double delta = T/N; //時間間隔

  double u = exp(sigma*sqrt(delta)); //上昇確率
  double d = exp(-sigma*sqrt(delta)); //下降確率
  double p = (exp(r*sqrt(delta)) - d) / (u - d); //リスク中立確率
  

/*株価ツリーの構築*/
  for(int i = 0; i < N+1; i++) {
    for(int j = 0; j < N+1; j++) {
      S[i][j] = S0*u^j*d^(i-j);
      
    }
  }

/*満期でのオプション価格*/
  for(int j = 0; j < N+1; j++) {
    if(op == 1) C[N][j] = max(S[N][j] - K, 0);
    if(op == 0) C[N][j] = max(K - S[N][j], 0);
    }

/*backwardにオプション価格を計算*/
  for(int i = N; i >= 0; i--) {
    for(int j = 0; j < N+1; j++) {
      C[i][j] = exp(-r*sqrt(delta))*(p*C[i+1][j+1] + (1-p)*C[i+1][j]);
    }
  }

  return C[0][0]  
}

int main() {
  
  //入力
  cin >> r >> sigma >> T >> K >> S0 >> N;
  print(binary_tree(0.1, 0.2, 5/12, 60, 62, 300, 1));
 // binary_tree(double r, double sigma, double T, double K, double S0, double N, int op = 1)
        
}                                                                                                           
