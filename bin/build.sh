#!/bin/bash

cd ..
cd build

# 运行 cmake 配置
cmake ..

# 编译项目
make

# 安装
make install
