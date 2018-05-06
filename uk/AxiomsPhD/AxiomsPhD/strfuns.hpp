#include <vector>
#include <string>

namespace strfuns {

    std::vector<std::string> split(const std::string& inStr, char sep);
    std::string trim(std::string& inStr); //Removes spaces from the beginning and the ending of a string.
    std::string int2str(int i);
}
