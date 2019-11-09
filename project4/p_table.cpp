#include <iostream>
#include <fstream>

using namespace std;

int pretty_tables(string file, char **str_of_int)
{
    fstream fs;
    fs.open (file, std::fstream::in);
    if(fs.is_open())
    {
        string row;
        int num_colms;
        int c = 0;
        int num_rows = 0;
        while(getline(fs, row))
        {
            size_t entry = row.find(",");
            while(entry != string::npos)
            {
                c++;
                if(c != 1)
                    cout << '|';
                cout << '\t' << row.substr(0, entry) << '\t';
                row = row.substr(entry + 1);
                entry = row.find(",");
            }
            cout << "|\t" << row << '\n';
            if(num_rows == 0)
                cout << '\n';
            num_rows++;
            num_colms = c;
        }
    } else
        cout << "unable to open file!";
    return 0;
}

int main(int argc, char **argv)
{
    // This main function is made to compile a testing program, in the Makefile, along with IVP
    if(argc != 2)
    {
        cout << "usage is p_test in.csv array_int (data points, should be sorted)\n";
        return 0;
    }
    return pretty_tables(argv[1], &argv[2]);
}