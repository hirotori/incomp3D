# Incompressible Fluid FLow Solver For Equil-Spaced Rectilinear Mesh
等間隔直交格子のための非圧縮性流体ソルバ.

## 現在の状態
移流項と粘性項の評価を選択可能にした. 陰解法には対応出来ていない.

## 概要
等間隔直交格子で流れ計算を行う. 

## 離散化手法
有限体積法(コロケート格子)
- 移流項: 1次精度風上法
- 粘性項: 2次精度中心差分法
- 時間精度: 
  - `solver_fs_t` : 1次精度Euler陽解法
  - `solver_fs_imp_t` : 粘性項のみクランク･ニコルソン法
  - その他, ユーザが各自で用意する.

## 解法
フラクショナルステップ法


## 解析手順
### ビルド方法
#### プログラムの構成
本ソルバはライブラリ群`lib`とメインプログラムに分けられている. `lib`にあるファイル群は数値シミュレーション
に必要なモジュール群であり, これらから静的ライブラリを生成する. 

メインプログラムは各々のケースに応じて用意し, コンパイル時に上記のライブラリとリンクすることで
実行ファイルを生成する.

#### ライブラリのビルド
`lib`のビルドにはcmakeを利用する. 

1. まず初めに, ビルド用のディレクトリを作成しておく. 

2. 次にconfigureを行う. コンパイラがgfortranの場合は以下のようにする. 

    `cmake .. -DCMAKE_Fortran_COMPILER=gfortran`

    コンパイラオプションを適用する場合は上に`-DCMAKE_BUILD_TYPE=`を指定する. デバッグモードは`debug`, リリースの場合は`release`とする.

3. cmakeによるconfigureが完了後, `make`コマンドを実行することでコンパイルを実施する. 終了後アーカイブ`libincomp3d.a`, `libincomp3d_util.a`が生成される.

#### メインプログラムのコンパイル
例題として3つのサンプルケースを用意している.ここでは例として`sample`のコンパイル方法を説明する.

`./src`にはメインプログラムが入っている. 中身は本ライブラリの使用例で, 解析ケース, アルゴリズム, シミュレーターの3つのユーザ派生型でメインプログラムを構成している. 

メインプログラムのビルドは`make`で実施することになっている. サンプルのMakefileを用意しているのでそちらを参考に.

**注意** ライブラリとメインプログラムのコンパイラは同じである必要がある. `lib`のコンパイルを`gfortran`で行った場合, メインプログラムのコンパイルも同じコンパイラで実施すること. そうしないとコンパイルエラーが出る.

### プログラムの実行

解析にあたって, まず計算条件設定ファイル`config.txt`を編集する.
#### メッシュの生成
ファイルの2~3行目.
```
33 33 2  !grid points (imax,jmax,kmax)
3.141592653 3.141592653 0.01  !length (x,y,z)     
```
`grid points` は計算領域の格子点の個数を指定する.

`length`は計算領域の寸法を指定する.
#### 計算条件
注釈通りに設定する. 
#### 境界条件
次節を参照.

### 境界条件
境界条件は圧力, 速度それぞれに指定する. 
現在指定可能な条件は境界値固定, 境界勾配固定, 周期境界, バッファ付き周期境界の4種類である. 

境界条件はプログラム内部で以下のように番号付けされている.
```
bc_fix_val = 1              !固定値境界
bc_fix_grad = 2             !固定勾配
bc_fix_periodic = 3         !周期境界
bc_fix_periodic_buffer = 4  !バッファ付き周期境界
```

境界条件はコンフィグファイル`config.txt`で指定する.
フォーマットは以下の通りとなる.
```
3 3 0.0 0.0 0.0 0.0  !B.C. for i = 1    face. [bc_u, bc_p, u, v, w, p]
3 3 0.0 0.0 0.0 0.0  !B.C. for i = imax face. [bc_u, bc_p, u, v, w, p]
3 3 0.0 0.0 0.0 0.0  !B.C. for j = 1    face. [bc_u, bc_p, u, v, w, p]
3 3 1.0 0.0 0.0 0.0  !B.C. for j = jmax face. [bc_u, bc_p, u, v, w, p]
2 2 0.0 0.0 0.0 0.0  !B.C. for k = 1    face. [bc_u, bc_p, u, v, w, p]
2 2 0.0 0.0 0.0 0.0  !B.C. for k = kmax face. [bc_u, bc_p, u, v, w, p]
```

**境界面の値(`u,v,w,p`)は, 固定値境界の場合は速度/圧力として, 固定勾配の場合は速度勾配/圧力勾配として扱われる**.

数値シミュレーションで扱う境界条件の組み合わせは限られている. 以下では主な組み合わせについて説明する.

### 1. 流入条件
境界面に固定された速度（一様流）を与える. 通常,圧力はNeumann境界となる.
#### 例: 一様流入
速度$\bm{u}=(1,0,0)$ で流入する条件の場合は以下のように指定する.
```
1 2 1.0 0.0 0.0 0.0 ![u_bc, p_bc, u,v,w,p]
```
圧力をゼロ勾配とするため, `p`に`0.0`を指定している.


### 2. 流出条件
対流流出条件を与える. ここでは流れ方向の速度の微分がゼロとする. 圧力は境界面で固定する.
```
2 1 0.0 0.0 0.0 1.0 ![u_bc, p_bc, u,v,w,p]
```

速度にゼロ勾配を与えるため, `u,v,w`を全てゼロにしている.

**注意1** 勾配規定条件は面に垂直な勾配として考える. そのため, 例えば`i=1,imax` の面に対しては境界速度成分は`u`のみが有効. ただ`v,w`にゼロ以外を入れると狂うので全てゼロにしてください.

**注意2** 圧力を任意に指定することが出来るが, 普通は圧力の初期値と対応させるべき.

### 3. 滑り無し条件
滑り無し条件を課す. 境界面で速度をゼロに固定する. 圧力は流入条件と同様にNeumann境界となる.
```
1 2 0.0 0.0 0.0 0.0 ![u_bc, p_bc, u,v,w,p]
```

**注意** 速度はゼロ以外でも可能. 物理的に適切な値を考えて入れてください.
### 4. 滑り条件

### 5. 周期境界条件
2つの境界面に指定することで周期境界条件を達成する. 

**注意** 可能なペアは i=1とi=imx, j=1とi=jmx, k=1とk=kmxの三通りのみ.

### 6. 非定常流出条件
まだ作ってない. 移流方程式に基づいて流出条件を設定する.

### 7. バッファ付き周期境界条件
周期境界として参照する**速度**を任意の場所に指定する.

**注意** 例えばi=1をこれにしたとき, i=imxについては周期境界ではなく流出条件とするべき. 例:

```
4 2 0.0 0.0 0.0 0.0  !B.C. for i = 1    face. [bc_u, bc_p, u, v, w, p]
2 1 0.0 0.0 0.0 1.0  !B.C. for i = imax face. [bc_u, bc_p, u, v, w, p]
```