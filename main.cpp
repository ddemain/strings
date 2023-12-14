#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// the alphabet is chars {32:127, 13, 10}; the delimiter is char(0)=NUL (std::strings are not NT);
constexpr char alphabet[] = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\r\n";
// chars are by default ascii-encoded 1-byte values;
constexpr int pk = 257;

//class d/d

struct Hit {
    int start;
    int length;
    float accuracy;
    Hit(int s, int l, float a = 1.0f) : start(s), length(l), accuracy(a) {};
};

class Match {
private:
    string base;
    string pattern;
    int indent;
    bool sorted;
    vector<Hit> hits;
public:
    Match(string b, string p, vector<Hit> h, int indent = 5, bool sorted = true);
    friend ostream& operator<<(ostream& os, Match& match);
};

ostream &operator<<(ostream &os, Match &m) {
    os << "string = \"" << m.base << "\";\n";
    os << "pattern = \"" << m.pattern << "\", " << m.hits.size() << " hits produced (" << ((m.sorted) ? "sorted" : "unsorted") << ") \n";
    for (const Hit& h : m.hits) {
        int pi = h.start - m.indent;
        int si = (int) m.base.size() - (h.start + h.length + m.indent);
        string pre = (pi > 0) ?
                     "..." + m.base.substr(pi, m.indent) :
                     m.base.substr(0, m.indent + pi);
        string suf = (si > 0) ?
                     m.base.substr(h.start + h.length, m.indent) + "..." :
                     m.base.substr(m.base.size() - (m.indent + si), m.indent + si);
        os << "hit (" << h.accuracy * 100 << "%, pos " << h.start << " to " << h.start + h.length << "): " << pre << "<" << m.base.substr(h.start, h.length) << ">" << suf << "\n";
    } //kind of done
    return os;
}

Match::Match(string b, string p, vector<Hit> h, int indent, bool sorted) : base(std::move(b)), pattern(std::move(p)) {
    if (sorted) sort(h.begin(),h.end(),[](const Hit& a, const Hit& b){
            return a.accuracy > b.accuracy; //hits sorted based on accuracy
        });
    this->sorted = sorted;
    this->indent = indent;
    hits = h;
}

//end of class d/d

Match Naive(const string& base, const string& pattern) { //computes in O((n-m+1)*m) ~ O(n^2)
    vector<Hit> hits;
    int n = (int) base.size();
    int m = (int) pattern.size();
    for (int i = 0; i < n-m; ++i) { // n-m times
        if (base.substr(i,m) == pattern) { // O(m), charwise
            hits.emplace_back(i, m);
        }
    }
    return {base, pattern, hits};
}

Match RabinKarp(const string& base, const string& pattern) { //computes in O(m+n-m+1) ~ O(n)
    vector<Hit> hits;
    int n = (int) base.size();
    int m = (int) pattern.size();
    long long int hp = 0, hb = 0, fmult = 1;
    for (int i = 0; i < m; ++i) { // O(m)
        hp = hp + pattern[i] * fmult; // += char*pk^i
        hb = hb + base[i] * fmult; // += char*pk^i
        fmult = fmult * pk; // pk^(i+1)
    } //fmult == pk^patternSize
    for (int i = 0; i < n - m; ++i) { // n-m times
        if (hb == hp) { hits.emplace_back(i, m); } // O(1)
        hb = (hb - base[i] + base[i + m] * fmult) / pk;
        // [strin]g = s*pk^0 + t*pk^1 + ... + n*pk^(m-1) ->
        // [0trin]g = 0 + t*pk^1 + ... + n*pk^(m-1) ->
        // [0tring] -> 0 + t*pk^1 + ... + n*pk^(m-1) + g*pk^m ->
        // [tring] -> t*pk^0 + ... n*pk(m-2) + g^pk(m-1)
    }
    return {base, pattern, hits};
}

//vector<int> prefixFunction (const string& s) {
//    int n = (int) s.size();
//    vector<int> prefixArray (n);
//    for (int i = 1; i < n; ++i) {
//        int j = prefixArray[i-1];
//        while (j > 0 && s[i] != s[j])
//            j = prefixArray[j-1];
//        if (s[i] == s[j]) ++j;
//        prefixArray[i] = j;
//    }
//    return prefixArray;
//}
//
//vector<vector<int>> transitionTable(string pattern) {
//    int n = (int) pattern.size();
//    vector<int> prefixArray(prefixFunction(pattern));
//    vector<vector<int>> table (n, vector<int> (98));
//    for (int i = 0; i < n; ++i) { //n ops
//        for (char c: alphabet) { //A ops
//            if (i > 0 && c != pattern[i])
//                table[i][c] = table[prefixArray[i-1]][c];
//            else
//                table[i][c] = i + (c == pattern[i]);
//        }
//    }
//    return table; // O(nA)
//}
//
//Match FiniteAutomaton(const string& base, const string& pattern) {
//    vector<Hit> res;
//    int baseSize = (int) base.size();
//    int patternSize = (int) pattern.size();
//    //transition computer
//}

Match KnuthMorrisPratt(const string& base, const string& pattern) {
    vector<Hit> hits;
    int m = (int) pattern.size();

    string s = pattern + char(0) + base; //composite string; NUL to restrict the prefix
    int l = (int) s.size();
    vector<int> prefix (l);

    for (int i = 1; i < l; ++i) { //n+m times
        int j = prefix[i-1];
        while (j > 0 && s[i] != s[j])
            j = prefix[j-1];
        if (s[i] == s[j]) ++j;
        prefix[i] = j;

        if (prefix[i] == m) {
            hits.emplace_back(i - 2 * m, m);
        }
    }
    return {base, pattern, hits};
}

Match BoyerMoore(const string& base, const string& pattern) { // gnu grep ftw
    vector<Hit> hits;
    int n = (int) base.size();
    int m = (int) pattern.size();
    for (int i = 0; i < n-m; ++i) { // n-m times
        if (base.substr(i,m) == pattern) { // O(m), charwise
            hits.emplace_back(i, m);
        }
    }
    return {base, pattern, hits};
}

int main() {
//    string x = "queLorem ipsum dolor sit amet, consectetur adipiscing elit. Morbi pellentesque rutrum mauris a pretium. Duis sodales vitae lorem id vulputate. Nullam vitae dui interdum, sollicitudin urna quis, mollis ligula. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Vestibulum turpis augue, cursus vel mi non, dictum convallis metus. Nunc et leo efficitur, auctor est in, porttitor libero. Ut vulputate cursus condimentum.\n"
//               "Vestibulum sit amet fermentum lorem, at dictum nunc. Aliquam scelerisque condimentum massa a blandit. Vestibulum eu velit sagittis, tincidunt dolor ac, iaculis lacus. Integer quis varius ligula. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Suspendisse ut fermentum libero, at pretium nisl. Pellentesque consectetur mi tortor, id elementum felis eleifend eu. Duis vehicula eget sapien eget ultrices. Ut sem lectus, pulvinar ac est sed, rutrum mattis purus. Duis ultricies enim accumsan ante finibus suscipit. Ut consectetur velit a eros commodo, sed iaculis neque vulputate. Nulla venenatis rhoncus porttitor. Pellentesque blandit venenatis felis, eleifend consequat mauris consectetur a. Praesent eget vulputate sapien. Sed rutrum cursus lectus id consequat.\n"
//               "In ac tortor at odio ornare posuere. Mauris gravida neque a diam sodales tempor. Quisque pellentesque lacus nisi, ac fermentum lacus rhoncus vel. Sed ac viverra orci. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus in convallis nulla. Nam faucibus nisi nec posuere pulvinar. Maecenas fringilla quam in ultricies scelerisque. Proin ac mi et ex malesuada dictum. Nullam tincidunt leo lacus, et porta sapien cursus porta. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus.\n"
//               "Mauris ultrices vel sem vel laoreet. Phasellus lacus nisl, tincidunt vel libero et, semper pharetra libero. Fusce augue diam, tristique ac blandit quis, finibus vel mi. Curabitur at dolor pretium ex ornare posuere. Mauris egestas eros neque, sit amet tempor sem ornare sed. Pellentesque turpis mi, tincidunt ut efficitur sed, feugiat sit amet magna. Maecenas quam lectus, iaculis nec elit non, vehicula convallis tortor. Aliquam viverra efficitur molestie. Fusce pulvinar ac odio sit amet dignissim. In vitae risus feugiat, convallis turpis a, ullamcorper leo. Donec ornare leo justo, ut interdum neque sodales sed. Donec elit nisi, congue eu metus quis, condimentum egestas orci. Nam eget sem nibh.\n"
//               "Nullam vitae enim ut odio ornare maximus. Aliquam malesuada felis ex, sed tristique diam egestas vel. Ut ac egestas elit, sed condimentum est. Suspendisse potenti. Fusce feugiat dictum lacus at lobortis. Sed consectetur nunc id pretium vestibulum. Morbi viverra mauris et dapibus mattis. Donec ultricies augue tincidunt mauris dignissim, ac commodo mauris suscipitque.";
//    string y = "que"
    string x = "Sampletextsamplestringsample.";
    string y = "ampl";
    Match naive = Naive(x, y);
    Match rk = RabinKarp(x, y);
//    Match fa = FiniteAutomaton(s, "ed");
    Match kmp = KnuthMorrisPratt(x, y);
    Match bm = BoyerMoore(x, y);
    cout << "Naive:\n" << naive << "\n";
    cout << "Rabin-Karp:\n" << rk << "\n";
    cout << "Knuth-Morris-Pratt:\n" << kmp << "\n";
    return 0;
}
