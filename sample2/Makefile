#==========================================================
# MakeFile for main application
# first version : 2022/09/19
# Author        : T.Ikeda
#==========================================================
# メインプログラム名.
SRC := ./src/fs3duns.f90

# ビルドされる実行プログラムの名称.
TARGET = fs3duns

#コンパイラ
FC := gfortran

#コンパイルオプション
FFLAGS := -O0 -fbacktrace
#----------------------------------------------------------
# 以下共通
#----------------------------------------------------------
# モジュールファイルの場所.
INCLUDE = ../build/mod
# コンパイル済みアーカイブファイルの場所.
LDFLAGS = -L../build/lib -lincomp3d -L../build/lib/common -lincomp3d_util

${TARGET}: ${SRC}
	${FC} ${SRC} -I${INCLUDE} ${LDFLAGS} -o $@ ${FFLAGS}

clean_output :
	rm output/*

clean :
	rm ${TARGET}