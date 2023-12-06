#pragma once

#include <fstream>
#include "../Geometry/Mesh.h"


class IO {

    public:
        IO();
        ~IO();

        static void WriteVTU(Mesh &mesh, std::string filename);



    private:
};
