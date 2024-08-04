#include "commonfunctions.h"

void putProcess(double procedure, double time_used) {

    std::cout << std::fixed << std::setprecision(3) << "Processing: " << procedure * 100 << "%" \
            << "\t\tTime: " << timeFormatting(time_used).str() \
            << "\t\tEstimate remaining time: " << timeFormatting((1 - procedure) / procedure * time_used).str() << std::endl;
    
}

std::stringstream timeFormatting(double microSeconds) {

    std::stringstream ret;
    ret << (unsigned long long)(microSeconds*1000.0) << "mius" << " (";
    unsigned long long seconds = microSeconds ;
    if (seconds < 60) {
        ret << seconds << "s";
    }
    else if (seconds < 3600) {
        ret << seconds / 60ull << "min " << seconds % 60ull << "s";
    }
    else {
        ret << seconds / 3600ull << "h " << seconds % 3600ull / 60ull << "min " << seconds % 60ull << "s";
    }
    ret << ")";
    return ret;

}

unsigned long long currentTime() {

    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::chrono::system_clock::duration duration = now.time_since_epoch();
    unsigned long long microSecondsOfDuration = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    return microSecondsOfDuration;

}

bool checksame_vector(vector <int> vec1, vector <int> vec2){

    sort(vec1.begin(), vec1.end());
    sort(vec2.begin(), vec2.end());
    if (vec1.size() != vec2.size()){
        return false;
    }
    for (int i = 0; i < vec1.size(); i++){
        if (vec1[i] != vec2[i]){
            return false;
        }
    }
    return true;

}

// template <typename tn> read(int &a){
//     tn x = 0, f = 1; char c = getchar();
//     for (; !isdigit(c); c = getchar()) if (c == '-') f = -1;
//     for (; isdigit(c); c = getchar()) x = x * 10 + c - 48;
//     a = x * f;
// }

void read(int &a){
    int x = 0, f = 1; char c = getchar();
    for (; !isdigit(c); c = getchar()) if (c == '-') f = -1;
    for (; isdigit(c); c = getchar()) x = x * 10 + c - 48;
    a = x * f;
}

void print_vec(vector <int> vec, string st, bool flag){
    if (flag == false) return;
    cout << st;
    for (int i = 0, si = vec.size(); i < si; i++)
        cout << vec[i] << (i == si - 1 ? '\n' : ' ');
}
