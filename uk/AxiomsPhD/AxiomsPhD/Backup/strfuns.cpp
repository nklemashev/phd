#include "strfuns.hpp"
#include <vector>
#include <string>

std::vector<std::string> strfuns::split(
    const std::string& inStr,
    char sep)
{
    std::vector<std::string> retVec(0);
    unsigned int i, j;
    char buf[1024];
    j = 0;
    for (i = 0; i < inStr.length(); i++)
    {
        if (inStr[i] == sep)
        {
            buf[j] = 0;
            std::string str(buf);
            retVec.push_back(str);
            j = 0;
        }
        else
        {
            buf[j++] = inStr[i];
        }
    }
    //Last string
    buf[j] = 0;
    std::string str(buf);
    retVec.push_back(str);
    return retVec;
}

std::string strfuns::trim(std::string& inStr)
{
    unsigned int i, j;
    //Beginning...
    for (i = 0; i < inStr.length(); i++)
    {
        if (inStr[i] == ' ')
        {
            continue;
        }
        else
        {
            break;
        }
    }
    // Now i indexes the first non-space character in the inStr
    if (i == inStr.length())
    {// inStr consists only of spaces
        return "";
    }
    for (j = inStr.length() - 1; j >= 0; j--)
    {
        if (inStr[j] == ' ')
        {
            continue;
        }
        else
        {
            break;
        }
    }
    // Now j indexes the last non-space character in the inStr
    // j cannot be < 0; j < 0 means that inStr consists only of spaces
    // and the function already returned empty string
    std::string retStr(inStr, i, j - i + 1);
    return retStr;
}

std::string strfuns::int2str(int i)
{
    std::string retStr;
    int ii, dig, iDiv;
    bool sign;
    if (i == 0)
    {
        retStr = "0";
        return retStr;
    }
    else if (i < 0)
    {
        retStr = "";
        ii = -1 * i;
        sign = false;
    }
    else //if i > 0
    {
        ii = i;
        retStr = "";
        sign = true;
    }

    while (ii > 0)
    {
        iDiv = ii / 10;
        dig = ii - iDiv * 10;
        retStr = (char)('0' + dig) + retStr;
        ii = iDiv;
    }

    if (sign)
    {
        return retStr;
    }
    else
    {
        return "-" + retStr;
    }
}
