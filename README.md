# ようこそMixtNormalSBパッケージへ
`MixtNormalSB`パッケージは混合正規分布まわりの計算、特に
* 乱数生成
* 確率密度関数の計算
* 対数尤度関数の計算
* EMアルゴリズムによる最尤推定およびクラスタリング

をサポートするR言語のパッケージです。主に混合正規分布の理解の助けになるパッケージを提供したいという思いから、株式会社すうがくぶんかの教材のひとつとして作成しています。なお、このパッケージではtidyverse styleを尊重するよう心掛けています。

## Introduction
`MixtNormalSB`パッケージは`devtools`パッケージの`install_github`関数を用いてインストールすることができます。
```
# Let's install and load our package !
install_github("utaka233/MixtNormalSB")
library(MixtNormalSB)
```
代表的な使い方の一つは、1次元混合正規分布からの乱数生成とEM algorithmによる最尤推定です。特に
* EM algorithmによって推定された母集団分布のパラメータ
* 各iterationごとの対数尤度の履歴
* 各データポイントのクラスタリングの結果

はすべて、最初から`tibble`の形でoutputされます。
```
# 乱数生成
library(dplyr)
library(ggplot2)
x <- random_mixt_normal(n = 100, mu = c(-4.0, 4.0), sigma = c(1.0, 4.0), pi = c(0.25, 0.75))
data <- data_frame(x = x)
ggplot(data = data, mapping = aes(x = x)) + geom_histogram(binwidth = 2.0)
```
```
# EM algorithm
result_em <- x %>%
    EM_mixt_normal(max_iter = 100,
                   tol = 1.0,
                   init_mu = c(1.0, 1.0),
                   init_sigma = c(1.0, 1.0),
                   init_pi = c(0.5, 0.5))
result_em$params    # 母集団パラメータの最尤推定値
result_em$estimated_component    # 各標本のクラスタリングの結果
result_em$log_likelihood_history    # 対数尤度の更新履歴
ggplot(data = result_em$log_likelihood_history,
       mapping = aes(x = iter, y = log_likelihood)) +
       geom_line() +
       ggtitle("History of log likelihood")
```
