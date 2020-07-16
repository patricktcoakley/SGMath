# Simple Game Math

## A simple, easy to use C++14 math library

The goal for Simple Game Math was to create a lightweight math library that could be dropped in as a way to do 2D or 3D game math without a lot of dependencies or files. As someone who is learning computer graphics and game development, I found some of the existing solutions to be somewhat unintuitive to use and had source code that wasn't easy to digest or debug; I wrote this as an exercise to both better understand some of the math involved in computer graphics and to have something that made sense to me as a user. 

To that end, the design objectives for Simple Game Math are to maintain its simplicity while gradually increasing performance and adding convenience features over time. This means that I will probably keep it non-templated, use minimal modern C++ features where practical, and maintain it as a single-header project. 

Simple Game Math is in active development, and there are probably (definitely) errors in some of the math. While I feel that I have a decent base version right now, the next immediate goal is to create tests and iterate over the existing interface a bit more.