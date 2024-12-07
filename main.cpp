#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;

//对于数量级为10
//No.1 快速
//No.2 堆
//No.3 归并

//对于数量级100
//No.1 基数
//No.2 快速
//No.3 归并

//更大的数位
//No.1 计数
//No.2 基数
//No.3 快速

//综合排名(10数量级权重0.05,100数量级权重0.15，高数量级权重0.8）
//No.1 基数( 3*0.05 + 7*0.15 + 6*0.8 = 6.00 )
//No.2 计数( 1*0.05 + 2*0.15 + 7*0.8 = 5.95 )
//No.3 快速( 7*0.05 + 6*0.15 + 5*0.8 = 5.25 )
//No.4 归并( 5*0.05 + 5*0.15 + 4*0.8 = 4.20 )
//No.5 堆  ( 6*0.05 + 4*0.15 + 3*0.8 = 3.30 )
//No.6 希尔( 4*0.05 + 3*0.15 + 2*0.8 = 2.25 )
//No.7 桶  ( 2*0.05 + 1*0.15 + 1*0.8 = 1.05 )


#define NUM 100000
#define MAX 100000

class SortResult {
public:
    string name;
    double time;
};

class Solution {
public:
    // 计数排序
    void sort1(vector<int>& arr, SortResult& result) {
        auto start = chrono::high_resolution_clock::now();
        int max = *max_element(arr.begin(), arr.end());
        vector<int> count(max + 1);
        for (auto i : arr) count[i]++;
        int index = 0;
        for (int i = 0; i <= max; i++)
            for (int j = 0; j < count[i]; j++)
                arr[index++] = i;
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "计数排序";
        result.time = elapsed.count();
    }

    // 堆排序
    void sort2(vector<int>& arr, SortResult& result) {
        auto start = chrono::high_resolution_clock::now();
        buildHeap(arr);
        for (int i = arr.size() - 1; i > 0; i--) {
            swap(arr[0], arr[i]);
            heapify(arr, i, 0);
        }
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "堆  排序";
        result.time = elapsed.count();
    }

    // 希尔排序
    void sort3(vector<int>& arr, SortResult& result) {
        auto start = chrono::high_resolution_clock::now();
        int n = arr.size();
        vector<int> gaps;
        int k = 0;
        while (true) {
            int gap = 9 * (1 << (2 * k++)) - 9 * (1 << k) + 1;
            if (gap >= n) break;
            gaps.push_back(gap);
        }
        reverse(gaps.begin(), gaps.end());
        for (int gap : gaps)
            for (int i = gap; i < n; i++) {
                int temp = arr[i], j = i;
                while (j >= gap && arr[j - gap] > temp) {
                    arr[j] = arr[j - gap];
                    j -= gap;
                }
                arr[j] = temp;
            }
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "希尔排序";
        result.time = elapsed.count();
    }

    // 基数排序
    void sort4(vector<int>& arr, SortResult& result) {
        auto start = chrono::high_resolution_clock::now();
        radixSort(arr);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "基数排序";
        result.time = elapsed.count();
    }

    // 归并排序
    void sort5(vector<int>& arr, SortResult& result) {
        auto start = chrono::high_resolution_clock::now();
        vector<int> temp(arr.size());
        mergeSortIterative(arr, temp, 0, arr.size() - 1);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "归并排序";
        result.time = elapsed.count();
    }

    // 快速排序
    void sort6(vector<int>& arr, SortResult& result) {
        auto start = chrono::high_resolution_clock::now();
        quickSort(arr, 0, arr.size() - 1);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "快速排序";
        result.time = elapsed.count();
    }

    // 桶排序
    void sort7(vector<int>& arr, SortResult& result) {
        if (arr.empty()) return;
        auto start = chrono::high_resolution_clock::now();

        int minElement = *min_element(arr.begin(), arr.end());
        int maxElement = *max_element(arr.begin(), arr.end());
        int range = maxElement - minElement + 1;

        vector<vector<int>> buckets(range);

        for (int num : arr) {
            buckets[num - minElement].push_back(num);
        }

        size_t index = 0;
        for (auto& bucket : buckets) {
            if (bucket.size() > 1) {
                insertionSort(bucket);
            }
            copy(bucket.begin(), bucket.end(), arr.begin() + index);
            index += bucket.size();
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        result.name = "桶  排序";
        result.time = elapsed.count();
    }

private:
    void heapify(vector<int>& arr, int n, int root) {
        int largest = root, left = 2 * root + 1, right = 2 * root + 2;
        if (left < n && arr[left] > arr[largest]) largest = left;
        if (right < n && arr[right] > arr[largest]) largest = right;
        if (largest != root) {
            swap(arr[root], arr[largest]);
            heapify(arr, n, largest);
        }
    }

    void buildHeap(vector<int>& arr) {
        int n = arr.size();
        for (int i = n / 2 - 1; i >= 0; i--) heapify(arr, n, i);
    }

    void insertionSort(vector<int>& arr, int left, int right) {
        for (int i = left + 1; i <= right; ++i) {
            int key = arr[i];
            int j = i - 1;
            while (j >= left && arr[j] > key) {
                arr[j + 1] = arr[j];
                --j;
            }
            arr[j + 1] = key;
        }
    }

    void merge(vector<int>& arr, vector<int>& temp, int left, int mid, int right) {
        int i = left, j = mid + 1, k = left;
        while (i <= mid && j <= right) {
            if (arr[i] <= arr[j]) {
                temp[k++] = arr[i++];
            }
            else {
                temp[k++] = arr[j++];
            }
        }
        while (i <= mid) {
            temp[k++] = arr[i++];
        }
        while (j <= right) {
            temp[k++] = arr[j++];
        }
        for (int p = left; p <= right; ++p) {
            arr[p] = temp[p];
        }
    }

    void mergeSortIterative(vector<int>& arr, vector<int>& temp, int left, int right) {
        int n = right - left + 1;
        for (int curr_size = 1; curr_size <= n - 1; curr_size *= 2) {
            for (int left_start = 0; left_start < n - 1; left_start += 2 * curr_size) {
                int mid = min(left_start + curr_size - 1, n - 1);
                int right_end = min(left_start + 2 * curr_size - 1, n - 1);
                if (mid < right_end) {
                    merge(arr, temp, left_start, mid, right_end);
                }
            }
        }
    }

    void quickSort(vector<int>& vec, int left, int right) {
        while (left < right) {
            int pivotIndex = medianOfThree(vec, left, right);
            int pivot = vec[pivotIndex];
            swap(vec[pivotIndex], vec[right]);
            int i = left - 1;
            for (int j = left; j < right; j++)
                if (vec[j] <= pivot) swap(vec[++i], vec[j]);
            swap(vec[i + 1], vec[right]);
            if (i - left < right - i - 1) {
                quickSort(vec, left, i);
                left = i + 2;
            }
            else {
                quickSort(vec, i + 2, right);
                right = i;
            }
        }
    }

    int medianOfThree(vector<int>& vec, int left, int right) {
        int mid = left + (right - left) / 2;
        if (vec[left] > vec[mid]) swap(vec[left], vec[mid]);
        if (vec[left] > vec[right]) swap(vec[left], vec[right]);
        if (vec[mid] > vec[right]) swap(vec[mid], vec[right]);
        return mid;
    }

    void countingSortByDigit(vector<int>& vec, vector<int>& output, int exp) {
        int n = vec.size();
        vector<int> count(10, 0);

        for (int i = 0; i < n; ++i)
            count[(vec[i] / exp) % 10]++;

        for (int i = 1; i < 10; ++i)
            count[i] += count[i - 1];

        for (int i = n - 1; i >= 0; --i) {
            int digit = (vec[i] / exp) % 10;
            output[count[digit] - 1] = vec[i];
            count[digit]--;
        }

        vec.assign(output.begin(), output.end());
    }

    void radixSort(vector<int>& vec) {
        if (vec.empty()) return;
        int maxElement = *max_element(vec.begin(), vec.end());
        vector<int> output(vec.size());

        for (int exp = 1; maxElement / exp > 0; exp *= 10)
            countingSortByDigit(vec, output, exp);
    }

    void insertionSort(vector<int>& arr) {
        for (size_t i = 1; i < arr.size(); ++i) {
            int key = arr[i];
            size_t j = i;
            while (j > 0 && arr[j - 1] > key) {
                arr[j] = arr[j - 1];
                --j;
            }
            arr[j] = key;
        }
    }
};

bool compare(SortResult a, SortResult b) {
    return a.time < b.time;
}

int main() {
    Solution s;
    vector<int> arr;
    std::random_device rd;  // 用于获取随机种子
    std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎
    std::uniform_int_distribution<> dis(1, MAX);
    for (int i = 0; i < NUM; i++) {
        arr.emplace_back(dis(gen));
    }
    vector<int> arr1(arr);
    vector<SortResult> results;
    results.push_back({ "计数排序", 0 });
    results.push_back({ "基数排序", 0 });
    results.push_back({ "快速排序", 0 });
    results.push_back({ "归并排序", 0 });
    results.push_back({ "堆  排序", 0 });
    results.push_back({ "希尔排序", 0 });
    results.push_back({ "桶  排序", 0 });

    s.sort1(arr1, results[0]);
    arr1 = arr;
    s.sort4(arr1, results[1]);
    arr1 = arr;
    s.sort6(arr1, results[2]);
    arr1 = arr;
    s.sort5(arr1, results[3]);
    arr1 = arr;
    s.sort2(arr1, results[4]);
    arr1 = arr;
    s.sort3(arr1, results[5]);
    arr1 = arr;
    s.sort7(arr1, results[6]);

    sort(results.begin(), results.end(), compare);

    cout << "七大排序竞赛现场，让我们拭目以待！" << endl;
    cout << "对" << NUM << "个1 ~ " << MAX << "的随机数进行排序" << endl;
    for (const auto& result : results) {
        cout << result.name << "耗时\t" << result.time << " s\n";
    }
    return 0;
}