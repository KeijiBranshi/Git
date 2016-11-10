#include <iostream>
#include <string>

bool isValidChar(char c) {
  bool isNum = (c>=48) && (c<=57);
  bool isCapitalLet = (c>=65) && (c<=90);
  bool isLowerLet = (c>=97) && (c<=122);

  return isNum || isCapitalLet || isLowerLet;
}
bool isPalindrome(std::string s) {
  int size = s.size();
  int i=0;
  int j=size-1;

  while (i<j) {
    if (!isValidChar(s[i])) {
      ++i; continue;
    }
    else if (!isValidChar(s[j])) {
      --j; continue;
    }
    if (tolower(s[i++]) != tolower(s[j--])) return false;
  }

  return true;
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "I need strrinnggsss" << std::endl;
    return -1;
  }
  else if (argc > 2) {
    std::cout << "Put it all in quotes pleeeassee" << std::endl;
    return -1;
  }

  std::string s(argv[1]);
  std::cout << std::boolalpha << isPalindrome(s) << std::endl;

  return 0;
}
