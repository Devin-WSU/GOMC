/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef WRITER_H
#define WRITER_H

#include <string> //for file names etc.
#include <stdlib.h> //for exit

#include <iostream>
#include <fstream>
#include <ostream>

class Writer
{
public:
  //Init later...
  Writer() {}
  //Init now.
  Writer(std::string const& name, std::string const& alias,
         const bool crit = true, const bool note = true,
         const bool isBin = false)
  {
    Init(name, alias, crit, note, isBin);
  }

  ~Writer(void)
  {
    if (isOpen) close();
  }

  //Set main class vars.
  void Init(std::string const& name, std::string const& alias,
            const bool crit, const bool note, const bool isBin = false)
  {
    fileName = name;
    fileAlias = alias;
    critical = crit;
    notify = note;
    isOpen = false;
    firstWrite = true;
    nameWAlias = fileAlias + " file: ./"  + fileName;
    isBinary = isBin;
  }

  //Open or close a file, with basic protections
  void open(void)
  {
    if (isOpen) return;
    file.open(fileName.c_str(),
              isBinary ? (std::ios::out | std::ios::binary) : (firstWrite ? std::ios::out : std::ios::app));
    CheckFileState(true, "...could not be opened.", "Writing to ");
  }


  //Open or close a file, with basic protections
  void openOverwrite(void)
  {
    if (isOpen) return;
    file.open(fileName.c_str(), isBinary ? (std::ios::out | std::ios::binary) : std::ios::out);
    CheckFileState(true, "...could not be opened.", "Writing to ");
  }


  void close(void)
  {
    if (!isOpen) return;
    file.close();
    CheckFileState(false, "...could not be closed.", "Finished writing ");
    //If first write is complete, unset firstWrite flag.
    if (firstWrite)
      firstWrite = false;
  }

  //Make public to expose object.
  std::ofstream file;
protected:

  bool GoodFileWData(void)
  {
    return file.is_open() && file.good();
  }

  void HandleError(std::string const& msg)
  {
    std::cerr << ((critical) ? "Error " : "Warning ") << nameWAlias << std::endl
              << msg << std::endl;
    if (critical)
      exit(1);
  }

  void HandleNote(std::string const& msg)
  {
    std::cout << msg << nameWAlias << std::endl;
  }

  void CheckFileState(const bool expected,
                      std::string const & errMessage, std::string const& note)
  {
    isOpen = GoodFileWData();
    if (isOpen == expected && notify)
      HandleNote(note);
    else if (isOpen != expected)
      HandleError(errMessage);
  }

  std::string fileName, fileAlias, nameWAlias;
  bool critical, notify, isOpen, firstWrite, isBinary;
};

#endif /*WRITER_H*/
