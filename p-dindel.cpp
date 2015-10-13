/**
 * p-dindel.cpp
 * This program calls indels using dindel in a multi-threaded fashion.
 * @author: Mohammad Shabbir Hasan
 * @email: shabbir5@cs.vt.edu
 * @version: 1.0
 * @Date: May 6, 2014
 * Virginia Tech, Blacksburg, USA. 
 */

#include <iostream>
#include <thread>
#include <math.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <time.h>

using namespace std;

double get_current_time(){
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

void *executeShell(void *threadId){
    long tid;
    tid = (long)threadId;

    stringstream shellFileName;
    shellFileName<<"tmp/thread_"<<tid<<".sh";

    stringstream permission_cmd;
    permission_cmd<<"chmod 777 "<<shellFileName.str()<<endl;
    system(permission_cmd.str().c_str());	

    stringstream execution_cmd;	
    execution_cmd<<"./"<<shellFileName.str()<<";";
    system((char*)execution_cmd.str().c_str());
}

void create_tmp_directory(){
    string cmd = "mkdir -p tmp";
    system(cmd.c_str());
}

int get_number_of_realignment_windows(){
    string cmd = "ls tmp | grep \"p-dindel_out_realign_windows\" > tmp/realignment_files_list.txt";
    
    try{
        system(cmd.c_str());
    }
    catch(...){
        cout<<"Exception in generating list of realignment windows"<<endl;
        return -1;
    }
    
    int number_of_realignment_windows = 0;
    string line;
    
    try{
        ifstream myfile("tmp/realignment_files_list.txt");

        while(getline(myfile, line)){
            number_of_realignment_windows++;
        }
        myfile.close();
    }
    catch(...){
        cout<<"Exception in reading the list of realignment windows"<<endl;
        return -1;
    }
    
    return number_of_realignment_windows;
}

int get_number_of_cores(){
    int number_of_cores = sysconf(_SC_NPROCESSORS_ONLN);
    
    return number_of_cores;
}

int get_number_of_threads_per_core(){
    int number_of_thread_per_core = 1;
    
    string cmd = "lscpu | grep \"Thread(s) per core:\" > tmp/thread_info.txt";
    system(cmd.c_str()); //system call to get the thread info
    
    string line, thread_info_line;
    
    try{
        ifstream thread_info_file ("tmp/thread_info.txt");
        if (thread_info_file.is_open()){
            while (getline(thread_info_file,line)){
                thread_info_line = line;
            }
        
            thread_info_file.close();
        }
                
        string arr[4];
        int i = 0;
        stringstream ssin(thread_info_line);
        while (ssin.good() && i < 4){
            ssin >> arr[i];
            i++;
        }
        
        number_of_thread_per_core = atoi(arr[3].c_str());
        
        if(number_of_thread_per_core <= 1){
            number_of_thread_per_core = 1;
        }
    }
    catch(...){
        cout<<"Exception in getting thread info"<<endl;
    }
    
    return number_of_thread_per_core;
}

void run_step1(string input, string reference){
    string step1_command = "./utils/dindel --analysis getCIGARindels --bamFile " + input + " --outputFile tmp/p-dindel_out --ref " + reference;
    try{
        system(step1_command.c_str());
    }
    catch(...){
        cout<<"Exception in running step 1"<<endl;
    }
}

void run_step2(){
    system("cp utils/makeWindows.py .");
    string step2_command = "./makeWindows.py --inputVarFile tmp/p-dindel_out.variants.txt --windowFilePrefix tmp/p-dindel_out_realign_windows --numWindowsPerFile 1000 ";
    try{
        system(step2_command.c_str());
    }
    catch(...){
        cout<<"Exception in running step 2"<<endl;    
    }
    system("rm makeWindows.py");
}

void run_step3(string input, string reference){
    int num_of_realignment_windows = get_number_of_realignment_windows();    
    int num_of_cores = get_number_of_cores();
    int num_of_threads_per_core = get_number_of_threads_per_core();
        
    int num_of_windows_per_thread = (int) ceil(num_of_realignment_windows/(double)(num_of_cores * num_of_threads_per_core));

    int windowIndex = 1;
    ofstream scriptFile[num_of_cores * num_of_threads_per_core];

    try{
        for(int i = 0; i<num_of_cores*num_of_threads_per_core; i++){
            stringstream fileName;
            fileName<<"tmp/thread_"<<i<<".sh";
            scriptFile[i].open(fileName.str());

            for(int j = windowIndex; j<=(windowIndex+num_of_windows_per_thread -1); j++){
                scriptFile[i]<<"./utils/dindel --analysis indels --doDiploid --bamFile " + input + " --ref "+ reference + " --varFile tmp/p-dindel_out_realign_windows."<<j<<".txt --libFile tmp/p-dindel_out.libraries.txt --outputFile tmp/p-dindel_out_realign_windows."<<j<<endl<<endl;			

                if(j>=num_of_realignment_windows)
                    break;
            }

            scriptFile[i].close();

            if(windowIndex >= num_of_realignment_windows)
                break;

            windowIndex+= num_of_windows_per_thread;		
        }
    }
    catch(...){
        cout<<"Exception in running step 3"<<endl;
    }
    
    pthread_t threads[num_of_cores*num_of_threads_per_core];
    int thread_creation_flag;

    for(int k = 0; k<num_of_cores*num_of_threads_per_core; k++){
        thread_creation_flag = pthread_create(&threads[k], NULL, executeShell, (void *)(intptr_t)k);

        if(thread_creation_flag){
            cout<<"Error: Unable to create thread, "<<thread_creation_flag<<endl;
            exit(-1);
        }
    }
    
    for(int l = 0; l<num_of_cores*num_of_threads_per_core; l++){
        pthread_join(threads[l], NULL);
    }    
}

void run_step4(){
    string step4_command = "ls tmp | sed -e 's/^/tmp\\//' |grep \".glf.txt\"  > p-dindel_out_indel_list.txt";
    
    try{
        system(step4_command.c_str());
    }
    catch(...){
        cout<<"Exception in running step 4"<<endl;
    }
}

void run_step5(string output, string reference){
    system("cp utils/mergeOutputDiploid.py .");
    string step5_command = "./mergeOutputDiploid.py --inputFiles p-dindel_out_indel_list.txt --outputFile "+output+"_p-dindel_output.vcf --ref "+reference;
    
    try{
        system(step5_command.c_str());
    }
    catch(...)
    {
        cout<<"Exception in running step 5"<<endl;
    }
    system("rm mergeOutputDiploid.py");
}

void run_clean_up(){
    try{
        system("rm -rf tmp");
        system("rm p-dindel_out_indel_list.txt");
        system("clear");
    }
    catch(...){
        cout<<"Exception in cleaning up step"<<endl;
    }            
}

int main(int argc, char* argv[]){    
    string input;
    string output;
    string reference;
    
    if(argc < 4){
        cout<<"Usage: "<<argv[0]<<" -i input.bam -o output_file_name -r reference.fa"<<endl;
        
        return 1;
    }
    
    else{        
        for(int n = 1; n<argc; n++){ // since at index 0, it is ./p-dindel i.e. name of the program.
            string arg = argv[n];
            
            if(arg == "-i"){
                input = argv[n+1];
            }
            
            if(arg == "-o"){
                output = argv[n+1];
            }
            
            if(arg == "-r"){
                reference = argv[n+1];
            }                        
        }
        
        cout<<endl<<endl<<"P-Dindel Starts Execution!!"<<endl<<endl;
        double start = get_current_time();
        create_tmp_directory();
        run_step1(input, reference);
        run_step2();
        run_step3(input, reference);
        run_step4();
        run_step5(output, reference);
        run_clean_up();
        double end = get_current_time();
        double execution_time = end - start;
        cout<<endl<<"P-Dindel Finishes Execution!!"<<endl<<endl;
        cout<<"Output saved at: "<<output<<"_p-dindel_output.vcf"<<endl<<endl;
        cout<<"Time Taken: "<<((int)(execution_time/60))<<" minutes "<<fmod(execution_time, 60)<<" seconds"<<endl<<endl;
        
        /****Will Delete it Later ****/
        string timing_file_name = output+"_p-dindel_timing.txt";
        ofstream timing_file;
        timing_file.open(timing_file_name);
        timing_file<<"Time Taken: "<<((int)(execution_time/60))<<" minutes "<<fmod(execution_time, 60)<<" seconds"<<endl<<endl;
        timing_file.close();
        /****Will Delete it Later ****/
    }    
    
    return 0;
}