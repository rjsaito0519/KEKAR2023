# (1)コンパイラ
CC=gcc
CXX=g++
# (2)コンパイルオプション
CFLAGS  = -O2

# (3)実行ファイル名
TARGET  = $(basename $(NAME))
# (4)コンパイル対象のソースコード
SRCS    = $(NAME)
# (5)オブジェクトファイル名
OBJDIR = ./obj/
OBJS    = $(patsubst %.C,$(OBJDIR)%.o,$(SRCS))

# (6)インクルードファイルのあるディレクトリパス
INCDIR  = -I./include

# (8)追加するライブラリファイル
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)


.PHONY: all
all: $(OBJDIR) $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)%.o: %.C
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $< $(CXXLIBS) $(CXXFLAGS) -lcurses -lSpectrum -ltbb

$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) $(INCLUDE) -o $@ $^ src/macro.cc src/nagao_macro.cc $(CXXLIBS) $(CXXFLAGS) -lcurses -lSpectrum -ltbb

# -lcursesは #include <ncurses.h> 用のライブラリ
# -lSpectrum は #include "TSpectrum.h" 用のライブラリ

.PHONY: clean
clean:
	rm -f *.d $(OBJDIR)*.o $(TARGET)
