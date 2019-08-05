#ifndef METAINFOLINE_H
#define METAINFOLINE_H
#include <string>
#include <vector>

class MetaInfoline
{
private:
  std::string fileformat;
  std::string fileDate;
  std::string source;
  std::string reference;

  std::vector<std::string> extra;

public:
  MetaInfoline();
  //get header
  std::string getDateStringNow();
  std::string getTextFileformat();
  std::string getTextFileDate();
  std::string getTextSource();
  std::string getTextReference();

  //set header
  void setheader_fileformat(std::string value);
  void setheader_fileDate(std::string value);
  void setheader_source(std::string value);
  void setheader_reference(std::string value);
  
};

#endif