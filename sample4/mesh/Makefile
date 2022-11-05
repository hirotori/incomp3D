#==========================================================
# MakeFile for main application
# first version : 2022/09/19
# Author        : T.Ikeda
#==========================================================
# メインプログラム名.
SRC := poiseuille_mesh.f90

# ビルドされる実行プログラムの名称.
TARGET = poiseuille_mesh

#コンパイラ
FC := gfortran

#コンパイルオプション
FFLAGS := -g -O0 -fbacktrace -fcheck=all

#モジュールファイルの出力場所
MODDIR = ./mod
#----------------------------------------------------------
# 以下共通
#----------------------------------------------------------
# モジュールファイルの場所.
INCLUDE = ../../build/mod
# コンパイル済みアーカイブファイルの場所.
LDFLAGS = -L../../build/lib -lincomp3d -L../../build/lib/common -lincomp3d_util

${TARGET}: ${SRC}
	${FC} ${SRC} -J ${MODDIR} -I${INCLUDE} ${LDFLAGS} -o $@ ${FFLAGS}

clean_output :
	rm output/*

clean :
	rm ${TARGET}