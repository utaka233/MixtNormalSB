# ようこそmixturesbパッケージへ
`mixturesb`パッケージは混合正規分布まわりの計算、特に
* `random_mixt_normal` : 乱数生成
* `density_mixt_normal` : 確率密度関数の計算
* `LL_mixt_normal` : 対数尤度関数の計算
* `EM_mixt_normal` : EMアルゴリズムによる最尤推定およびクラスタリング 

をサポートするR言語のパッケージです。主に混合正規分布の理解の助けになるパッケージを提供したいという思いから、**株式会社すうがくぶんか**の教材のひとつとして作成しています。なお、このパッケージではtidyverse styleを尊重するよう心掛けています。また、将来的には混合分布全般の乱数生成やEM-algorithm、Bayes推定の提供を目標にしています。

## Installation
`mixturesb`パッケージは、`devtools`パッケージの`install_github`関数を用いてインストールすることができます。
```
# Let's install and load our package !
devtools::install_github("utaka233/mixturesb")
library(mixturesb)
```
代表的な使い方の一つは、1次元混合正規分布からの乱数生成とEM algorithmによる最尤推定です。特に
* EM algorithmによって推定された母集団分布のパラメータ
* 各iterationごとの対数尤度の履歴
* 各データポイントのクラスタリングの結果

をすべて`tibble`の形で取得することも、簡単に`summary`関数を用いて概要をチェックすることも可能です。`plot`関数を用いれば、EM algorithmの対数尤度の更新履歴を簡単に確認することもできます。

## 1. Usage
### 1.1 generate random values of 1-dim. mixtures of normal distribution
1次元混合正規分布からの乱数生成は次のように行うことが出来ます。
```
# 乱数生成
library(dplyr)
library(ggplot2)
x <- random_mixt_normal(n = 100, mu = c(-4.0, 4.0), sigma = c(1.0, 4.0), ratio = c(0.25, 0.75))
data <- tibble(x = x)
ggplot(data = data, mapping = aes(x = x)) + geom_histogram(binwidth = 2.0)
```
<img src="https://github.com/utaka233/garage/blob/master/imgs_mixturesb/histogram_of_x.png" alt = "ヒストグラム" width="250" />

### 1.2 EM algorithm for 1-dim. mixtures of normal distribution
先ほどサンプリングした1次元混合正規分布の乱数に対して、(log-)EM algorithmを用いて最尤推定を行ってみましょう。
```
# EM algorithm
result_em <- x %>%
    EM_mixt_normal(max_iter = 100,
                   tol = 0.01,
                   init_mu = c(-1.0, 1.0),
                   init_sigma = c(1.0, 1.0),
                   init_ratio = c(0.5, 0.5))
summary(result_em)    # 推定結果の概要
plot_LL(result_em)    # ggplot2によるiterationごとの対数尤度の更新履歴
plot_components(result_em, binwidth = 2.0)    # ggplot2によるcomponent推定の結果の可視化
```
<img src="https://github.com/utaka233/garage/blob/master/imgs_mixturesb/history_LL.png" alt = "EMアルゴリズムの更新履歴" width="200" /><img src="https://github.com/utaka233/garage/blob/master/imgs_mixturesb/plot_components.png" alt = "EMアルゴリズムの更新履歴" width="500" />

## 2. Getting tibble of the EM-algorithm results
```
# 各推定結果の詳細
result_em$params    # 母集団パラメータの最尤推定値に関するtibble
result_em$estimated_component    # 各標本のクラスタリングの結果のtibble
result_em$log_likelihood_history    # 対数尤度の更新履歴のtibble
```
