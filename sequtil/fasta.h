# include       <cstdio>
# include       <string>

/// Thanks to Peter "sehr Ocker" for this great quick hack!
class fasta_reader
{

protected:

	bool file_open;
	
public:
	
	/// open for reading from stdin
	fasta_reader ();
		
	fasta_reader (const std::string &filename);

	void open(const std::string &filename);
	
	~fasta_reader ();
	
	int get_seq (std::string &header, std::string &seq); // 0 = EOF, 1 = seq ok

	static void reverse_complement(std::string& seq);
	
private:
	
    enum        State { Void, Header, Data };
    FILE        *fp_;
    State       state_;
	 
};

