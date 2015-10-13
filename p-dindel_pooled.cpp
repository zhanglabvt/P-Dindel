/**
 * p-dindel_pooled.cpp
 * This program calls indels using dindel for pooled samples in a multi-threaded fashion.
 * @author: Mohammad Shabbir Hasan
 * @email: shabbir5@cs.vt.edu
 * @version: 1.0
 * @Date: May 30, 2014
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
    shellFileName<<"tmp_pooled/thread_"<<tid<<".sh";

    stringstream permission_cmd;
    permission_cmd<<"chmod 777 "<<shellFileName.str()<<endl;
    system(permission_cmd.str().c_str());	

    stringstream execution_cmd;	
    execution_cmd<<"./"<<shellFileName.str()<<";";
    system((char*)execution_cmd.str().c_str());
}

void create_tmp_directory(){
    string cmd = "mkdir -p tmp_pooled";
    system(cmd.c_str());
}

int get_number_of_realignment_windows(){
    string cmd = "ls tmp_pooled | grep \"p-dindel_pooled_out_realign_windows\" > tmp_pooled/realignment_files_list.txt";
    
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
        ifstream myfile("tmp_pooled/realignment_files_list.txt");

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
    
    string cmd = "lscpu | grep \"Thread(s) per core:\" > tmp_pooled/thread_info.txt";
    system(cmd.c_str()); //system call to get the thread info
    
    string line, thread_info_line;
    
    try{
        ifstream thread_info_file ("tmp_pooled/thread_info.txt");
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

void run_step1(string input, int sample_index, string reference){
    stringstream step1_command;
    step1_command<<"./utils/dindel --analysis getCIGARindels --bamFile "<<input<<" --outputFile tmp_pooled/"<<sample_index<<"_p-dindel_pooled_out --ref "<<reference;
    try{
        system(step1_command.str().c_str());
    }
    catch(...){
        cout<<"Exception in running step 1" <<endl;
    }
}

void run_step2(){
    system("cp utils/makeWindows.py .");
    string merge_all_samples_variant_file_command = "cat tmp_pooled/*.variants.txt > tmp_pooled/all_variants.txt";
    try{
        system(merge_all_samples_variant_file_command.c_str());
    }
    catch(...){
        cout<<"Exception in merging all samples variant files"<<endl;
    }
    
    string step2_command = "./makeWindows.py --inputVarFile tmp_pooled/all_variants.txt --windowFilePrefix tmp_pooled/p-dindel_pooled_out_realign_windows --numWindowsPerFile 1000 ";
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

    stringstream copy_bam_files_list_command;
    copy_bam_files_list_command<<"cp "<<input<<" tmp_pooled/all_bam_files.txt";
    try{
        system(copy_bam_files_list_command.str().c_str());
    }
    catch(...){
        cout<<"Exception in copying the list of all bam files"<<endl;
    }
    string merge_all_library_file_command = "cat tmp_pooled/*.libraries.txt > tmp_pooled/all_library.txt";
    try{
        system(merge_all_library_file_command.c_str());
    }
    catch(...){
        cout<<"Exception in merging library files of all samples"<<endl;
    }
    
    int windowIndex = 1;
    ofstream scriptFile[num_of_cores * num_of_threads_per_core];

    try{
        for(int i = 0; i<num_of_cores*num_of_threads_per_core; i++){
            stringstream fileName;
            fileName<<"tmp_pooled/thread_"<<i<<".sh";
            scriptFile[i].open(fileName.str());

            for(int j = windowIndex; j<=(windowIndex+num_of_windows_per_thread -1); j++){
                scriptFile[i]<<"./utils/dindel --analysis indels --doPooled --bamFiles tmp_pooled/all_bam_files.txt --ref "<<reference<<" --varFile tmp_pooled/p-dindel_pooled_out_realign_windows."<<j<<".txt --libFile tmp_pooled/all_library.txt --outputFile tmp_pooled/p-dindel_pooled_out_realign_windows."<<j<<endl<<endl;			

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
    string step4_command = "ls tmp_pooled | sed -e 's/^/tmp_pooled\\//' | grep \".glf.txt\"  > p-dindel_pooled_out_indel_list.txt";
    
    try{
        system(step4_command.c_str());
    }
    catch(...){
        cout<<"Exception in running step 4"<<endl;
    }
}

void run_step5(string output, string reference, int num_of_samples, int num_of_bam_files){
    system("cp utils/mergeOutputPooled.py .");
    stringstream step5_command;
    step5_command<<"./mergeOutputPooled.py --inputFiles p-dindel_pooled_out_indel_list.txt --outputFile "<<output<<"_p-dindel_pooled_output.vcf --ref "<<reference<<" --numSamples "<<num_of_samples<<" --numBamFiles "<<num_of_bam_files;
    try{
        system(step5_command.str().c_str());
    }
    catch(...)
    {
        cout<<"Exception in running step 5"<<endl;
    }
    system("rm mergeOutputPooled.py");
}

void run_clean_up(){
    try{
        system("rm -rf tmp_pooled");
        system("rm p-dindel_pooled_out_indel_list.txt");
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
    int number_of_bam_files = 0;
    int number_of_samples = 0;    
    
    if(argc < 4){
        cout<<"Usage: "<<argv[0]<<" -i input_bam_files_list -o output_file_name -r reference.fa"<<endl;
        
        return 1;
    }
    
    else{        
        for(int n = 1; n<argc; n++){ // since at index 0, it is ./p-dindel_pooled i.e. name of the program.
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
        
        cout<<endl<<endl<<"P-Dindel Pooled Starts Execution!!"<<endl<<endl;
        double start = get_current_time();
        create_tmp_directory();
        string bam_file_array[1000];
        int bam_file_index = 0;
        string line;
        try{
            ifstream myfile(input);

            while(getline(myfile, line)){
                bam_file_array[bam_file_index] = line;
                bam_file_index++;
            }
            myfile.close();
        }
        catch(...){
        cout<<"Exception in reading the list of bam files"<<endl;        
        }
        number_of_bam_files = bam_file_index - 1;
        number_of_samples = number_of_bam_files;        
        for(int i = 0; i<number_of_bam_files; i++){
            run_step1(bam_file_array[i], i, reference);        
        }        
        run_step2();
        run_step3(input, reference);
        run_step4();
        run_step5(output, reference, number_of_bam_files, number_of_samples);
        run_clean_up();
        double end = get_current_time();
        double execution_time = end - start;
        cout<<endl<<"P-Dindel Pooled Finishes Execution!!"<<endl<<endl;
        cout<<"Output saved at: "<<output<<"_p-dindel_pooled_output.vcf"<<endl<<endl;
        cout<<"Time Taken: "<<((int)(execution_time/60))<<" minutes "<<fmod(execution_time, 60)<<" seconds"<<endl<<endl;
        
        /****Will Delete it Later ****/
        string timing_file_name = output+"_p-dindel_pooled_timing.txt";
        ofstream timing_file;
        timing_file.open(timing_file_name);
        timing_file<<"Time Taken: "<<((int)(execution_time/60))<<" minutes "<<fmod(execution_time, 60)<<" seconds"<<endl<<endl;
        timing_file.close();
        /****Will Delete it Later ****/
    }    
    
    return 0;
}