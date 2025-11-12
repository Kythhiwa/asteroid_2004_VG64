#include <iostream>
#include "include/reader_obs.h"

int main()
{
   auto a = read_obs("data/ObsCodes.html");
   for (auto now : a)
   {
       std::cout << now.first << " " << now.second.Long << " " << now.second.cos << " " << now.second.sin << "\n";
   }
}
