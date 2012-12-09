
/// Thanks to Peter "sehr Ocker" for this great quick hack!


#include <iostream>
#include <algorithm>
#include "fasta.h"
#include "bbq_util.h"

// ==============================================================
fasta_reader::~fasta_reader()
// ==============================================================
 {
	 if (file_open)
		 fclose (fp_);
	 file_open = false;
 }

// ==============================================================
fasta_reader::fasta_reader()
// ==============================================================
    :   fp_ (0), state_ (Void)
{
	fp_ = stdin;
	file_open = false;
	if (!fp_)
		throw bbq_exception("Error reading from stdin.");
}   

// ==============================================================
fasta_reader::fasta_reader (const std::string &filename)
// ==============================================================
    :   fp_ (0), state_ (Void)
{
    fp_ = fopen (filename.c_str(), "rb");
	 file_open = true;
	 if (!fp_)
		 throw bbq_exception("Cannot open fasta file " + filename + " for reading. Try using pipes.");
}   

// ==============================================================
void fasta_reader::open (const std::string &filename)
// ==============================================================
{
	//std::cerr<<"OPEN FASTA\n";
	std::cerr.flush();
	if (file_open)
	{
		//std::cerr<<"FILE OPEN\n";
		std::cerr.flush();
		fclose (fp_);
	}
	
    fp_ = fopen (filename.c_str(), "rb");
	 if (!fp_)
		 throw bbq_exception("Cannot open fasta file " + filename + " for reading. Try using pipes.");
	 file_open = true;
}   



// ==============================================================
int fasta_reader::get_seq (std::string &header, std::string &seq)
// ==============================================================	
{
    int         result = -1;    // 0 = eof, 1 = good seq.

    //if (state_ != Void)
    //    MFATAL ("state_ != VOID");

    header.clear ();
    seq.clear ();

    while (result == -1)
    {
        int     c = fgetc (fp_);
		  if (!fp_)
			  throw bbq_exception("Error reading fasta file. Try using pipes.");

        switch (state_)
        {
            case Void:
                switch (c)
                {
                    case EOF:
                        result = 0;
                        state_ = Void;
                        break;

                    case ' ': case '\t': case '\n': case '\r':
                        break;

                    case '>':
                        state_ = Header;
                        header.push_back (c);
                        break;
            
                    default:
							  throw bbq_exception("illegal character in fasta file.\n");
                        break;
                }
                break;

            case Header:
                switch (c)
                {
                    case EOF:
                        result = 0;
                        break;

                    case '\n': case '\r':
                        state_ = Data;
                        break;

                    case '\001':    // ^A or Ctrl-A
							  throw bbq_exception("header contains Ctrl-A character (after %s)");
                        break;
    
                    default:
                        header.push_back (c);
                        break;
                }
                break;

            case Data:
					if (c=='>') {
						state_ = Void;
						ungetc (c, fp_);
						if (!fp_)
							throw bbq_exception("Error reading fasta file. Try using pipes.");
						result = 1;
					} else if (c==EOF) {
						result = 1;
						state_ = Void;
					} else if ((c>='a' && c<='z') || (c>='A' && c<='Z')) {
						seq.push_back(c);
					}
                break;

            default:
					throw bbq_exception("error reading fasta file.\n");
					break;
        }
    }

    return result;
}


// ==============================================================
void fasta_reader::reverse_complement(std::string& seq)
// ==============================================================	
{
	reverse(seq.begin(), seq.end());
	
	std::string::iterator s_it = seq.begin();
	std::string::iterator s_end = seq.end();
	for (; s_it!=s_end; ++s_it)
	{
		switch (*s_it)
		{
			case 'a': {
				*s_it = 't';
				break;
			}
			case 'A': {
				*s_it = 'T';
				break;
			}
			case 'c': {
				*s_it = 'g';
				break;
			}
			case 'C': {
				*s_it = 'G';
				break;
			}
			case 'g': {
				*s_it = 'c';
				break;
			}
			case 'G': {
				*s_it = 'C';
				break;
			}
			case 't': {
				*s_it = 'a';
				break;
			}
			case 'T': {
				*s_it = 'A';
				break;
			}
			case 'u': {
				*s_it = 'a';
				break;
			}
			case 'U': {
				*s_it = 'A';
				break;
			}
		}
	}
}
