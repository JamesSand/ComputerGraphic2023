#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

#include "render.hpp"
#include "scene.hpp"
using namespace std;
int main(int argc, char* argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum]
                  << std::endl;
    }
    if (argc < 4) {
        std::cout << "Usage: ./bin/PA1 <input scene file> <output bmp file> "
                     "<numRounds> <numPhotons> <ckpt_interval>"
                  << endl;
        return 1;
    }

    // 根据输入构建场景    
    Scene scene(argv[1]);

    cout << "done scenece parse" << endl;

    int numRounds = atoi(argv[3]);
    int numPhotons = atoi(argv[4]);
    int ckpt_interval = atoi(argv[5]);

    // 构建 sppm，arg2 是 output path
    SPPM sppm(scene, numRounds, numPhotons, ckpt_interval, argv[2]);
    sppm.render();


    return 0;
}
