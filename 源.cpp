#include<iostream>
#include<vector>
#include<stack>
#include<queue>
#include<set>
#include<map>
#include<numeric>
#include<string>
#include<iterator>
#include<list>
#include<algorithm>
#include<unordered_map>
#include<deque>
#include<unordered_map>
#include<unordered_set>
using namespace std;

struct TreeNode {
	int val;
	struct TreeNode *left;
	struct TreeNode *right;
	TreeNode(int x) :
		val(x), left(NULL), right(NULL) {
	}
};
struct B {
	TreeNode * p;
	int row;
	B(TreeNode*pp, int r) :p(pp), row(r) {}
};
vector<vector<int> > Print(TreeNode* pRoot) {
	queue<B> q;
	vector<vector<int> > r;
	if (!pRoot)return r;

	q.push(B(pRoot, 0));
	while (!q.empty()) {
		B t = q.front();
		r.resize(t.row + 1);
		if (t.row % 2 == 0)r[t.row].push_back(t.p->val);
		else r[t.row].insert(r[t.row].begin(), t.p->val);
		if (t.p->left)q.push(B(t.p->left, t.row + 1));
		if (t.p->right)q.push(B(t.p->right, t.row + 1));
		q.pop();
	}
	return r;
}
vector<int> twoSum(vector<int>& nums, int target) {
	multimap<int, int> m;
	for (int i = 0; i < nums.size(); i++)
		m.insert(pair<int, int>(nums[i], i));

	//auto p = m.equal_range(9);
	for (int i = 0; i < nums.size(); i++) {
		vector<int> r;
		auto p = m.equal_range(target - nums[i]);
		for (auto j = p.first; j != p.second; j++) {
			if (i != j->second) {
				r.push_back(i);
				r.push_back(j->second);
				return r;
			}
		}
	}
}
int expandAroundCenter(string s, int left, int right) {
	while (left >= 0 && right < s.length()
		&& s.at(left) == s.at(right)){
		left--;
		right++;
	}
	return right - left - 1;
}
string longestPalindrome(string s) {
	if (s.length() <= 1)return s;
	int st = 0, e = 0;
 	int l = 0;
	for (int i = 0; i < s.length(); i++) {
		int l1 = expandAroundCenter(s,i, i);
		int l2 = expandAroundCenter(s,i, i + 1);
		int len = max(l1, l2);
		if (len > e - st) {
			st = i - (len - 1) / 2;
			e = i + len / 2;
		}
	}
	return	s.substr(st, e-st+1);
}
double findMedianSortedArrays(vector<int> A, vector<int>B) {
	int m = A.size();
	int n = B.size();
	if (m > n) { // to ensure m<=n
		vector<int> temp = A; A = B; B = temp;
		int tmp = m; m = n; n = tmp;
	}
	int iMin = 0, iMax = m, halfLen = (m + n + 1) / 2;
	while (iMin <= iMax) {
		int i = (iMin + iMax) / 2;
		int j = halfLen - i;
		if (i<iMax && B[j - 1] > A[i]) {
			iMin = i + 1; // i is too small
		}
		else if ( i>iMin && A[i - 1] > B[j]) {
			iMax = i - 1; // i is too big
		}
		else { // i is perfect
			int maxLeft = 0;
			if (i == 0) { maxLeft = B[j - 1]; }
			else if (j == 0) { maxLeft = A[i - 1]; }
			else { maxLeft = max(A[i - 1], B[j - 1]); }
			if ((m + n) % 2 == 1) { return maxLeft; }

			int minRight = 0;
			if (i == m) { minRight = B[j]; }
			else if (j == n) { minRight = A[i]; }
			else { minRight = min(B[j], A[i]); }

			return (maxLeft + minRight) / 2.0;
		}
	}
	return 0.0;
}
string convert(string s, int numRows) {
	vector<string> strs(numRows);
	int z = 0;
	for (int i = 0; i < s.size(); i++){
		if (z >= 2 * numRows - 2){
			z = 0;
			strs[0] += s[i];}
		else if (z >= numRows)strs[2*numRows-2-z]+=s[i];
		else strs[z] += s[i];
		z++;
	}
	string r;
	for (int i = 0; i < numRows; i++) {
		r += strs[i];
	}
	return r;
}
int reverse(int x) {
	stack<int>s;
	int r = 0;
	while (x != 0) {
		s.push(x % 10);
		x = x / 10;
	}
	while (!s.empty()) {
		if (r > INT_MAX / 10 || r == INT_MAX / 10 && s.top() > 7)return 0;
		if (r < INT_MIN / 10 || r == INT_MIN / 10 && s.top() < -8)return 0;
		r = r * 10 + s.top();
		s.pop();
	}
	return r;
}
int myAtoi(string str) {
	if (str.size()<= 0)return 0;
	int i = 0;
	while (str[i] == ' ')i++;
	int sign = 1;
	if (str[i] == '-') { sign = -1; i++; }
	else if (str[i] == '+')i++;
	int r = 0;
	while (i != str.length()) {
		char c = str.at(i);
		if (c <= '9' && c >= '0') {
			if (r > INT_MAX / 10 || r == INT_MAX / 10 && c>'7')return INT_MAX;
			if (r < INT_MIN / 10 || r == INT_MIN / 10 && c >'8')return INT_MIN;
			r = r * 10 + sign*(int)(c - 48);
			i++;
		}
		else 
			return r;
	}
	return r;
}
bool isPalindrome(int xx) {
	int x = xx;
	if (x<0)return false;
	queue<int>s;
	int r = 0;
	while (x != 0) {
		s.push(x % 10);
		x = x / 10;
	}
	while (!s.empty()) {
		if (r > INT_MAX / 10 || r == INT_MAX / 10 && s.front() > 7)return 0;
		if (r < INT_MIN / 10 || r == INT_MIN / 10 && s.front() < -8)return 0;
		r = r * 10 + s.front();
		s.pop();
	}
	return r == x;
}
bool isMatch(string s, string p) {
	if (s.size()!=0 && p.size()==0)return false;
	if (s.size() == 0 && p.size() == 0)return true;
	if (p.size() >= 2 && p.at(1) == '*') {
		if (s.size() != 0)
			if (s.at(0) == p.at(0) || p.at(0) == '.')
				return isMatch(s.substr(1, s.length() - 1),
					p) ||
				isMatch(s, p.substr(2, p.length() - 2));
			else
				return isMatch(s, p.substr(2,p.length() - 2));
		else
			return isMatch(s, p.substr(2,p.length() - 2));
	}
	else {
		if (s.size() != 0)
			if (s.at(0) == p.at(0) || p.at(0) == '.')
				return isMatch(s.substr(1, s.length() - 1),
					p.substr(1, p.length() - 1));
			else
				return false;
		else
			return false;
	}
	
}
string intToRoman(int num) {
	vector<vector<string>> t = { {"","I","II","III","IV","V","VI","VII","VIII","IX",}, 
								 { "","X","XX","XXX","XL","L","LX","LXX","LXXX","XC" } ,
								 {"","C","CC","CCC","CD","D","DC","DCC","DCCC","CM"},
								 {"","M","MM","MMM"}};
	//int t = 0;
	int i = 0;
	string r;
	while (num != 0) {
		r = t[i][num % 10]+r;
		i++;
		num /= 10;
	}
	return r;
}
int romanToInt(string s) {
	int r = 0;
	int i = 0;
	while (i < s.size()) {
		if (s[i] == 'V') { r += 5; i++; }
		if (s[i] == 'L') {
			r += 50; i++;}
		if (s[i] == 'D') {
			r += 500; i++;}
		if (s[i] == 'I')
			if (i < s.size() - 1 && s[i + 1] == 'X'){
				r += 9;
				i += 2;	}
			else if (i < s.size() - 1 && s[i + 1] == 'V'){
				r += 4;
				i += 2;}
			else{
				r += 1;
				i += 1;}
		if (s[i] == 'X')
			if (i < s.size() - 1 && s[i + 1] == 'C') {
				r += 90;
				i += 2;
			}
			else if (i < s.size() - 1 && s[i + 1] == 'L') {
				r += 40;
				i += 2;
			}
			else {
				r += 10;
				i += 1;
			}
		if (s[i] == 'C')
			if (i < s.size() - 1 && s[i + 1] == 'M') {
				r += 900;
				i += 2;
			}
			else if (i < s.size() - 1 && s[i + 1] == 'D') {
				r += 400;
				i += 2;
			}
			else {
				r += 100;
				i += 1;
			}
		if (s[i] == 'M')
		{
			r += 1000;
			i += 1;
		}
	}
	return r;
}

string longestCommonPrefix(vector<string>& strs) {
	string lcp;
	int i = 0;
	if (strs.size() == 0)return lcp;
	while (1) {
		char c;
		if (i >= strs[0].size())return lcp;
		else c = strs[0][i];
		int j = 1;
		for (; j < strs.size(); j++) {
			if (i>=strs[i].size() || 
				strs[j][i] != c)return lcp;
		}
		if (j == strs.size())lcp += c;
		i++;
	}
}
vector<vector<int>> threeSum(vector<int>& nums) {
	vector<vector<int>> rr;
	for (int i = 0; i < nums.size(); i++) {
		map<int, int> m;
		for (int j = i+1; j < nums.size(); j++) {
			//if (j == i)continue;
			vector<int> r;
			if (m.find(-nums[i] - nums[j]) != m.end()) {
				r.push_back(nums[i]);
				r.push_back(nums[j]);
				r.push_back(-nums[i] - nums[j]);
				sort(r.begin(), r.end());
				if(find(rr.begin(),rr.end(),r)==rr.end())rr.push_back(r);
				m[nums[j]] = j;
			}
			else {
				m[nums[j]] = j;
			}
		}}
	return rr;
}
vector<vector<int>> threeSum2(vector<int> nums) {
	vector<vector<int>> res;
	sort(nums.begin(),nums.end());
	for (int i = 0; i + 2 < nums.size(); i++) {
		if (i > 0 && nums[i] == nums[i - 1]) {              // skip same result
			continue;
		}
		int j = i + 1, k = nums.size() - 1;
		int target = -nums[i];
		while (j < k) {
			if (nums[j] + nums[k] == target) {
				vector<int> t; t.push_back(nums[i]); t.push_back(nums[j]); t.push_back(nums[k]);
				res.push_back(t);
				j++;
				k--;
				while (j < k && nums[j] == nums[j - 1]) j++;  // skip same result
				while (j < k && nums[k] == nums[k + 1]) k--;  // skip same result
			}
			else if (nums[j] + nums[k] > target) {
				k--;
			}
			else {
				j++;
			}
		}
	}
	return res;
}
int threeSumCloset(vector<int> nums,int target) {
	//vector<vector<int>> res;
	int closet;
	closet = nums[0] + nums[1] + nums[2];
	sort(nums.begin(), nums.end());
	for (int i = 0; i + 2 < nums.size(); i++) {//i=0~n-3
		int j = i + 1, k = nums.size() - 1;
		int tar = target-nums[i];
		//closet = nums[i] + nums[j] + nums[k];
		while (j < k) {
			if (abs(nums[i] + nums[j] + nums[k] - target) < abs(closet - target))
				closet = nums[i] + nums[j] + nums[k];
			/*while (j<k -1&& nums[j+1] == nums[j])j++;
			while (j<k-1 && nums[k - 1] == nums[k])k--;*/
			if (nums[j] + nums[k] == tar) 
				return target;
			else if (nums[j] + nums[k] > tar) 
				k--;
			else 
				j++;
			
		}
	}
	return closet;
}
vector<string> letterCombinations(string digits) {
	vector<string> m = { "abc","def","ghi","jkl","mno","pqrs","tuv","wxyz" };
	vector<int> f(digits.size(), 0);
	vector<string> r;
	while (1) {
		string t;
		for (int i = 0; i < digits.size(); i++)t += m[digits[i] - 48 - 2][f[i]];
		r.push_back(t);
		int up = 0;
		f[0]+=1;
		for (int i = 0; i < digits.size(); i++) {
			f[i] += up;
			up = 0;
			if (f[i] >= m[digits[i] - 48 - 2].size()) {
				f[i] = 0; up = 1;}
			if (i == digits.size() - 1 && up == 1)goto finish;
		}
	}
finish:
	return r;
}
void dfs(string r, int level, int in, int out,int n,vector<string>& ret) {
	if (level == 2 * n) {
		ret.push_back(r); return;}
	if (in < n)dfs(r+'(', level+1, in + 1, out, n,ret);
	if (in > 0 && in >out)dfs(r + ')', level + 1, in , out+1, n,ret);
}
vector<string> generateParenthesis(int n) {
	vector<string> ret;
	dfs("(", 1, 1, 0, n, ret);
	return ret;
}
int longestValidParentheses(string s) {
	stack<int> stack;
	vector<pair<int, int>> valid;
	for (int i = 0; i < s.size(); i++) {
		if (s[i] == '(')stack.push(i);
		if (s[i] == ')') {
			if (!stack.empty()) {
				valid.push_back(pair<int,int>(stack.top(), i));
				while (1) {
					int j = valid.size() - 1;
					if (j < 1)break;
					if (valid[j - 1].first>valid[j].first) {
						swap(valid[j - 1], valid[j]);
						valid.pop_back();
					}else if (valid[j - 1].second+1==valid[j].first) {
						valid[j - 1].second = valid[j].second;
						valid.pop_back();					}
					else break;
				}
				stack.pop();
			}
			else{
				while (!stack.empty())stack.pop();
			}}}
	int max = 0;
	for (int i = 0; i < valid.size(); i++) {
		int t = valid[i].second - valid[i].first + 1;
		if (t > max)max = t;
	}
	return max;
}
int longestValidParentheses2(string s) {
	int left = 0, right = 0, maxlength = 0;
	for (int i = 0; i < s.length(); i++) {
		if (s[i] == '(') {
			left++;
		}
		else {
			right++;
		}
		if (left == right) {
			maxlength =max(maxlength, 2 * right);
		}
		else if (right >= left) {
			left = right = 0;
		}
	}
	left = right = 0;
	for (int i = s.length() - 1; i >= 0; i--) {
		if (s[i] == '(') {
			left++;
		}
		else {
			right++;
		}
		if (left == right) {
			maxlength = max(maxlength, 2 * left);
		}
		else if (left >= right) {
			left = right = 0;
		}
	}
	return maxlength;
}
struct ListNode {
	int val;
	ListNode *next;
	ListNode(int x) : val(x), next(NULL) {}
};
ListNode* reverseKGroup(ListNode* head, int k) {
	if (!head)return NULL;
	ListNode* dum = new ListNode(-1); dum->next = head;
	ListNode* p, *pp, *ppp, *h = dum;
	while(1){
		p = h->next;
		if (!p)return dum->next;//下一段一个都没了
		pp = p->next;
		if(!pp)return dum->next;//下一段只剩一个
		//下一段有两个及以上
		ppp = pp->next;
		//如果不足k个，直接返回了
		ListNode* tt = p;
		for (int i = 0; i < k; i++) {
			tt = tt->next;
			if (!tt)return dum->next;
		}
		for (int i = 1; i <= k - 1; i++) {
			pp->next = p;
			p = pp;
			pp = ppp;
			if (ppp)ppp = ppp->next; 
			else {//当前段已经没有了
				h->next->next = NULL;//本段第一个（将变为转置后最后一个）的next指向NULL
				h->next = p;//将上一段最后一个指向本段原最后一个
				return dum->next;
			}}
		h->next->next = pp;//本段第一个（将变为转置后最后一个）的next指向下一段第一个
		ppp = h->next;//临时存储本段第一个
		h->next = p;//上一段最后一个指向本段最后一个
		h = ppp;//将指示上段最后一个的变量设为指向本段最后一个（原本段第一个）
	}
}
ListNode* MakeList2(vector<int> nums) {
	if (nums.size() == 0)return NULL;
	ListNode* h=new ListNode(nums[0]),*p=h;
	for (int i = 1; i < nums.size(); i++) {
		p->next = new ListNode(nums[i]);
		p = p->next;
	}
	return h;
}
int removeDuplicates(vector<int>& nums) {
	if (nums.size() <= 1)return nums.size();
	int i = 0;
	while (i<nums.size() - 1) {
		if (nums[i + 1] == nums[i])nums.erase(nums.begin() + i + 1);
		else i++;
	}
	return nums.size();
}
vector<int> findSubstring(string s, vector<string>& words) {
	vector<int> permu(words.size());
	set<int> r;
	vector<int>rr;
	for (int i = 0; i < permu.size(); i++)permu[i] = i;
	do {
		string t;
		for (int i = 0; i < permu.size(); i++)t = t + words[permu[i]];
		int pos = -1;
		do {
			pos = s.find(t,pos+1);
			if (pos != string::npos)r.insert(pos); 
			else break;
		} while (1);
	} while (next_permutation(permu.begin(), permu.end()));
	for (auto i = r.begin(); i != r.end(); i++)rr.push_back(*i);
	return rr;
}
vector<vector<int>> combinationSum(vector<int>& candidates, int target) {
	vector<vector<int>> r;
	if (candidates.empty())return r;
	vector<int> stack;
	//stack.push_back(0);
	
	//bool pop = false;
	sort(candidates.begin(), candidates.end());
	stack.push_back(candidates[0]);
	int sum = candidates[0];
	while (1) {
		if(sum>=target){
			if (sum == target) {
				vector<int>out(stack.begin(), stack.end());
				r.push_back(out);
				//goto pop;
			}
			vector<int>::iterator index;
			while(1){
				sum -= stack.back();
				stack.pop_back();
				if (stack.empty()) goto combinationSum_end;
				int top = stack.back();
				index = find(candidates.begin(), candidates.end(), top);
				//若本次top仍有可以替换的candidates[i]，则替换，否则将上一层的元素也pop掉
				if (index + 1 != candidates.end())break;
			}
			sum -= stack.back();
			stack.back() = *(index + 1);
			sum += stack.back();
		}
		else {
			stack.push_back(stack.back());
			sum += stack.back();
		}
	}
combinationSum_end:
	return r;
}
vector<vector<int>> combinationSum2(vector<int>& candidates, int target) {
	vector<vector<int>> r;
	struct stack_frame {
		int num;
		int level;
		stack_frame(int n,int l):num(n),level(l){}
	};
	if (candidates.empty())return r;
	vector<stack_frame> stack;
	vector<int> res;
	sort(candidates.begin(), candidates.end());
	for (auto i = candidates.rbegin(); i != candidates.rend(); i++){
		stack.push_back(stack_frame(*i,1));
	}
	while (!stack.empty()) {
		//获取当前处理节点
		stack_frame t = stack.back();
		stack.pop_back();
		//计算和，以决定是
		while (res.size() >= t.level)res.pop_back();
		res.push_back(t.num);
		int sum = accumulate(res.begin(), res.end(), 0);
		if (sum == target)
			r.push_back(res); 
		else {
			//剪枝并加入后续节点
			int ub = upper_bound(candidates.begin(), candidates.end(), target - sum)-candidates.begin();
			int lb = lower_bound(candidates.begin(), candidates.end(), t.num)-candidates.begin();
			if (ub - lb<=0) 
				continue;
			for (int i = ub - 1; i >= lb; i--) 
				stack.push_back(stack_frame(candidates[i],res.size()+1));
		}
	}
	return r;
}
vector<vector<int>> fourSum(vector<int>& nums, int target) {
	vector<vector<int>> r;
	if (nums.size()<4)return r;
	sort(nums.begin(), nums.end());
	for (int i = 0; i <= nums.size() - 4; i++){
		if (i > 0 && nums[i - 1] == nums[i])
			continue;
		for (int j = i + 1; j <= nums.size() - 3; j++) {
			if (j>i+1 && nums[j] == nums[j - 1])
				continue;	
			int tar = target-nums[i]-nums[j];
			int k = j + 1;
			int l = nums.size() - 1;
			while(k<l){
				if (nums[k] + nums[l]>tar)
					l--;
				else if (nums[k] + nums[l] < tar)
					k++;
				else {
					vector<int>vtmp;
					vtmp.push_back(nums[i]); vtmp.push_back(nums[j]); vtmp.push_back(nums[k]); vtmp.push_back(nums[l]);
					r.push_back(vtmp);
					k++;
					l--;
					while (k<l && nums[k] == nums[k - 1]) k++;
					while (k<l && nums[l] == nums[l + 1]) l--;
					//break;
				}}}}
	return r;
}
void combinationSum2_dfs(vector<int>& candidates, vector<int> res, int index, vector<vector<int>>&r, const int& target) {
	if (target == 0) {
		r.push_back(res);
		return;
	}
	//if (target < 0)return;
	for (int i = index; i < candidates.size(); i++) {
		if (i > index && candidates[i] == candidates[i - 1])
			continue;
		if (candidates[i]>target)return;
		res.push_back(candidates[i]);
		combinationSum2_dfs(candidates, res, i + 1, r, target - candidates[i]);
		res.pop_back();
	}
}
vector<vector<int>> combinationSumII(vector<int>& candidates, int target) {
	vector<vector<int>> r;
	if (candidates.empty())return r;
	vector<int>res;
	sort(candidates.begin(), candidates.end());
	combinationSum2_dfs(candidates, res,0, r, target);
	return r;
}
bool LeetCode_44_isMatch(string s, string p) {
	//下标0表示空串，i表示0到i-1的字串是否match
	vector<vector<bool>> dp(s.size()+1, vector<bool>(p.size()+1, false));
	
	dp[0][0] = true;//both s and p is empty means matched
	for (int i = 1; i <= p.size(); i++)
		if (p[i-1] == '*')dp[0][i] = dp[0][i - 1];//if the string is empty and last pattern is '*'
	
	for (int i = 1; i <= s.size(); i++)//注意p[]s[]的下标
		for (int j = 1; j <= p.size(); j++) {
			if (s[i-1] == p[j-1] || p[j-1] == '?')dp[i][j] = dp[i - 1][j - 1];
			//* matches only this character or nothing or this character and maybe next
			else if (p[j-1] == '*')dp[i][j] = dp[i - 1][j - 1] || dp[i][j - 1]|| dp[i-1][j];
		}
	return dp[s.size()][p.size() ];
}
vector<vector<string>> LeetCode49_groupAnagrams(vector<string>& strs) {
	//用一个26位的数组记录每个单词中对应字母出现的次数来作map的key
	vector<vector<string>> ret;
	//ret.reserve(100);
	if (strs.empty())return ret;
	int key[26];
	map<string, int> map;
	for (string str : strs) {
		//初始化key
		memset(key, 0, sizeof(key));//sizeof数组即为数组元素*元素大小
		//计算key
		for (char c : str) 
			key[c - 'a']++;
		string str_key;
		for (int i = 0; i <= 25; i++)
			str_key = str_key + "#" + to_string(key[i]);
		//若key存在和不存在分别处理
		if (map.find(str_key) != map.end()) {
			int index = map[str_key];
			ret[index].push_back(str);
		}
		else {
			map[str_key] = ret.size();
			vector<string> tmp(1,str);
			//tmp->push_back(str);
			ret.push_back(tmp);
		}
	}
	return ret;
}
bool place_ith_queen(vector<int>& place, int i) {
	//给定的摆法，找寻第i个queen的下一个位置，若找到返回true，否则返回false
	int n = place.size();
	if (place[i] == n - 1)return false;//摆到最后一格了，没地方摆了
	for (int j = place[i]+1; j < n; j++) {
		//尝试将queen放在i行j列
			//检测前面放的queen是否在同一列或同一对角线
			int k = 0;
			for (; k < i; k++)
				if (place[k] == j || k-place[k] ==i - j || place[k] + k == j+i)
					break;//注意判断是否在同一斜线的条件
			//此摆法可行，返回true
			if(k==i){
				place[i] = j; return true;}
	}
	//没有下一个摆法了
	return false;
}
vector<vector<string>> LeetCode52_solveNQueens(int n) {
	vector<int> place(n, -1);
	vector<vector<string>> ret;
	//处理n<=4的情况
	if (n == 1) {
		vector<string> tmp;
		tmp.push_back("Q");
		ret.push_back(tmp);
		return ret;
	}
	if (n<4)return ret;
	int i = 0;
	bool b;
	//当不是第一个queen都摆不下的情况（第一个queen摆到最右了）
	while (!(!(b=place_ith_queen(place, i)) && i==0)) {
		if (b){
			i++;//摆放成功，摆下一个
			//最后一个摆放成功了
			if (i == n) {
				//输出解
				vector<string> t_vs;
				//构造一个解方阵，t代表第a行
				for (int a = 0; a < n; a++) {
					string t(n, '.');
					t[place[a]] = 'Q';
					t_vs.push_back(t);
				}
				ret.push_back(t_vs);
				place[--i] = -1;//当前行解法重置0
				i--;//重摆上一个
			}
		}
		//摆放不成功，上一个重摆
		else {
			place[i] = -1;//当前行解法重置0
			i--;
		}
	}
	return ret;
}
int maxSubArray(vector<int>& nums) {
	vector<int>dp(nums.size(), 0);
	dp[0] = nums[0];
	int max = dp[0];
	for (int i = 1; i < nums.size(); i++) {
		if (dp[i - 1] >= 0)dp[i] = dp[i - 1] + nums[i];
		else dp[i] = nums[i];
		if (dp[i] > max)max = dp[i];
	}
	return max;
}
vector<int> printMatrix(vector<vector<int> > matrix) {
	vector<int> o;
	if (matrix.empty())return o;
	int i = 0, j = 0;
	int u = -1, d = matrix.size(), l = -1, r = matrix[0].size();//四个方向的界
	while (1) {
		//往右走
		for (; j<r; j++)
			o.push_back(matrix[i][j]);
		j--;//将列标减回到行最后一个，并将行标下移
		i++;
		u = u + 1;
		if (o.size() >= matrix.size()*matrix[0].size())break;//所有元素已遍历完，走完了
		for (; i<d; i++)
			o.push_back(matrix[i][j]);
		i--;
		j--;
		r = r - 1;
		if (o.size() >= matrix.size()*matrix[0].size())break;
		for (; j>l; j--)
			o.push_back(matrix[i][j]);
		j++;
		i--;
		d = d - 1;
		if (o.size() >= matrix.size()*matrix[0].size())break;
		for (; i>u; i--)
			o.push_back(matrix[i][j]);
		i++;
		j++;
		l = l + 1;
		if (o.size() >= matrix.size()*matrix[0].size())break;
	}
	return o;
}
bool LeetCode55_canJump(vector<int>& nums) {
	int maxreach = 0;//上一步所能到达最远距离
	//int cur = 0;//dangqian
	int count = 0;//已走步数
	int i = 0;
	while (1) {
		count++;
		int t = maxreach;
		while (i <= t) {
			maxreach = max(maxreach, nums[i] + i);
			i++;
		}
		if (maxreach >= nums.size() - 1)return true;
		if (maxreach <= t)return false;
	}
}
struct Interval {
	int start;
	int end;
	Interval() : start(0), end(0) {}
	Interval(int s, int e) : start(s), end(e) {}
};
bool merge_comp(Interval intvl1, Interval intvl2) {
	return intvl1.start < intvl2.start || intvl1.start == intvl2.start && intvl1.end < intvl2.end;
}
vector<Interval> LeetCode56_merge(vector<Interval>& intervals) {
	sort(intervals.begin(), intervals.end(), merge_comp);
	int i = 0;
	int t;//此处intervals.size()被编译器优化，必须强制每次重算
	while (i <= (t=intervals.size()-2)) {
		if (intervals[i + 1].start <= intervals[i].end) {
			intervals[i].end = max(intervals[i + 1].end,intervals[i].end);
			intervals.erase(intervals.begin() + i + 1);
		}
		else
			i++;
	}
	return intervals;
}
vector<Interval> LeetCode57_insert(vector<Interval>& intervals, Interval newInterval) {
	int i = 0, j = intervals.size() - 1;
	//二分查找插入位置
	while (i <= j) {
		int mid = (i + j) / 2;
		if (intervals[mid].start > newInterval.start)
			j = mid - 1;
		else if (intervals[mid].start < newInterval.start)
			i = mid + 1;
		else {
			i = mid;
			break;
		}
	}
	intervals.insert(intervals.begin() + i, newInterval);
	i = 0;
	int t;//此处intervals.size()被编译器优化，必须强制每次重算
	while (i <= (t = intervals.size() - 2)) {
		if (intervals[i + 1].start <= intervals[i].end) {
			intervals[i].end = max(intervals[i + 1].end, intervals[i].end);
			intervals.erase(intervals.begin() + i + 1);
		}
		else
			i++;
	}
	return intervals;
}
vector<vector<int>> LeetCode59_generateMatrix(int n) {
	vector<vector<int>> matrix(n,vector<int>(n,0));
	if (n == 0)return matrix;
	int i = 0, j = 0;
	int u = -1, d = n, l = -1, r = n;//四个方向的界
	int count = 0;//当前数字
	while (1) {
		//往右走
		if (count >= n*n)break;//所有元素已遍历完，走完了
		for (; j < r; j++)
			matrix[i][j] = ++count;
		j--;//将列标减回到行最后一个，并将行标下移
		i++;
		u = u + 1;
		
		if (count >= n*n)break;//所有元素已遍历完，走完了
		for (; i<d; i++)
			matrix[i][j] = ++count;
		i--;
		j--;
		r = r - 1;
		
		if (count >= n*n)break;//所有元素已遍历完，走完了
		for (; j>l; j--)
			matrix[i][j] = ++count;
		j++;
		i--;
		d = d - 1;
		
		if (count >= n*n)break;//所有元素已遍历完，走完了
		for (; i>u; i--)
			matrix[i][j] = ++count;
		i++;
		j++;
		l = l + 1;
	}
	return matrix;
}
string LeetCode60_getPermutation(int n, int k) {
	//n!向量和数字向量
	vector<int> nding(n + 1, 1); for (int i = 1; i <= n; i++)nding[i] = nding[i - 1] * i;
	vector<int> nums(n); for (int i = 0; i < n; i++)nums[i] = i + 1;
	string ret;
	int i = n - 1;
	k = k - 1;
	while (!nums.empty()) {
		//从高位开始，每次从nums抽取对应的数字，直到nums为空
		int index = k / nding[i];
		ret += to_string(nums[index]);
		nums.erase(nums.begin() + index);
		k = k%nding[i];
		i--;
	}
	return ret;
}
ListNode* LeetCode61_rotateRight(ListNode* head, int k) {
	if (!head)return NULL;
	int len = 0;
	ListNode* p = head,*tail=NULL;
	while (p) {
		if (!p->next)tail = p;
		p = p->next;
		len++;
	}
	tail->next = head;
	p = head;
	k = len-k%len;//转换成左移次数
	//现在用tail保存上一个指针
	for (int i = 1; i <= k; i++) {
		tail = tail ->next;
		p = p->next;
	}
	tail->next = NULL;
	return p;
}
bool IsNum(char ch)
{
	if (ch<'0' || ch>'9') return false;
	else return true;
}

bool isNumeric(char* string)
{
	int i = 0;
	if (string[i] == '+' || string[i] == '-' || IsNum(string[i])) {
		while (string[++i] != '\0' && IsNum(string[i]));
		if (string[i] == '.') {
			if (IsNum(string[++i])) {
				while (string[++i] != '\0' && IsNum(string[i]));
				if (string[i] == 'e' || string[i] == 'E') {
					i++;
					if (string[i] == '+' || string[i] == '-' || IsNum(string[i])) {
						while (string[++i] != '\0' && IsNum(string[i]));
						if (string[i] == '\0') return true;
						else return false;
					}
					else return false;
				}
				else if (string[i] == '\0') return true;
				else return false;
			}
			else if (string[++i] == '\0') return true;
			else return false;
		}
		else if (string[i] == 'e' || string[i] == 'E') {
			i++;
			if (string[i] == '+' || string[i] == '-' || IsNum(string[i])) {
				while (string[++i] != '\0' && IsNum(string[i]));
				if (string[i] == '\0') return true;
				else return false;
			}
			else return false;
		}
		else if (string[i] == '\0') return true;
		else return false;
	}
	else return false;
}


enum state { Sta, Term, AM1, Mi1, N1, Dot1, Dot2, N2, E, AM2, N3 };
bool isNumeric2(char* string)
{
	state s = Sta;
	int j = 0;
	int i = 0;
	while (s != Term) {
		switch (s) {
		case Sta:
			if (string[i] == '-' || string[i] == '+')
				s = AM1;
			else if (string[i] == '.')
				s = Dot2;
			else if (isdigit(string[i]))
				s = N1;
			else return false;
			break;
		case AM1:
			if (isdigit(string[i]))s = N1;
			else if (string[i] == '.')s = Dot2;
			else return false;
			break;
		case N1:
			if (isdigit(string[i])) s = N1;
			else if (string[i] == '.') s = Dot1;
			else if (string[i] == '\0') s = Term;
			else if (string[i] == 'e' || string[i] == 'E') s = E;
			else return false;
			break;
		case Dot1:
			if (isdigit(string[i])) s = N2;
			else if (string[i] == '\0')s = Term;
			else if (string[i] == 'e' || string[i] == 'E') s = E;
			else return false;
			break;
		case Dot2:
			if (isdigit(string[i])) s = N2;
			else return false;
			break;
			break;
		case N2:
			if (isdigit(string[i])) s = N2;
			else if (string[i] == 'e' || string[i] == 'E') s = E;
			else if (string[i] == '\0') s = Term;
			else return false;
			break;
		case E:
			if (string[i] == '+' || string[i] == '-')
				s = AM2;
			else if (isdigit(string[i])) s = N3;
			else return false;
			break;
		case AM2:
			if (isdigit(string[i]))
				s = N3;
			else return false;
			break;
		case N3:
			if (isdigit(string[i])) s = N3;
			else if (string[i] == '\0') s = Term;
			else  return false;
			break;
		}
		i++;
	}
	return true;
}
bool match(char* str, char* pattern)
{
	if (*str == '\0' && *pattern == '\0')
		return true;
	if (*str != '\0' && *pattern == '\0')
		return false;
	//if the next character in pattern is not '*'
	if (*(pattern + 1) != '*')
	{
		if (*str == *pattern || (*str != '\0' && *pattern == '.'))
			return match(str + 1, pattern + 1);
		else
			return false;
	}
	//if the next character is '*'
	else
	{
		if (*str == *pattern || (*str != '\0' && *pattern == '.'))
			return match(str, pattern + 2) || match(str + 1, pattern);
		else
			return match(str, pattern + 2);
	}
}
bool match2(char* str, char* pattern) {
	//if (*str == '\0' && *pattern != '\0')return false;
	if (*str != '\0' && *pattern == '\0')return false;
	if (*str == '\0' && *pattern == '\0')return true;
	if (*(pattern + 1) == '*') {
		if (*str == *pattern || (*str != '\0' && *pattern == '.'))
			return match2(str + 1, pattern) || match2(str, pattern + 2);
		else
			return  match2(str, pattern + 2);
	}
	else {
		if (*str == *pattern || (*str != '\0' && *pattern == '.'))return match2(str + 1, pattern + 1);
		else return false;
	}
}
//struct ListNode {
//	int val;
//	ListNode *next;
//	ListNode(int x) : val(x), next(NULL) {}
//
//};
ListNode* swapPairs(ListNode* head) {
	if (!head)return NULL;
	ListNode* p = head;
	ListNode* pp = p->next;
	ListNode* ppp;
	if (pp) {
		head = pp;
		p->next = pp->next;
		pp->next = p;
		ppp = p;
		p = p->next;
		if (p)pp = p->next;
		else return head;
	}
	else return p;

	while (pp) {
		ListNode* t = pp->next;
		ppp->next = pp;
		pp->next = p;
		p->next = t;

		ppp = p;
		p = p->next;
		if (p)pp = p->next;
		else return head;
	}
	return head;
}
ListNode* MakeList(vector<int> arr) {
	if (arr.size() == 0)return NULL;
	ListNode* h = new ListNode(arr[0]);
	ListNode* p = h;
	for (int i = 1; i < arr.size(); i++) {
		p->next = new ListNode(arr[i]);
		p = p->next;
	}
	return h;
}
void nextPermutation(vector<int>& nums) {
	if (nums.size() <= 1)return;
	int partition = -1;
	for (int i = nums.size() - 2; i > 0; i--) {
		if (nums[i] < nums[i + 1]) {
			partition = i;
			break;
		}
	}
	if (partition == -1) { sort(nums.begin(), nums.end()); return; }
	for (int i = nums.size() - 1; i >= 0; i--) {
		if (nums[i] > nums[partition]) {
			swap(nums[i], nums[partition]);
			reverse(nums.begin() + partition + 1, nums.end());
			break;
		}
	}
	return;
}
int findfirst(vector<int>& nums, int start, int end, int target) {
	if (start > end)return -1;
	int mid = (start + end) / 2;
	//if (start == end && nums[mid] == target)return start; else return -1;
	if (nums[mid] > target)return findfirst(nums, start, mid - 1, target);
	else if (nums[mid] < target)return findfirst(nums, mid + 1, end, target);
	else if (mid - 1 >= 0 && nums[mid - 1] == target)return findfirst(nums, start, mid - 1, target);
	else return mid;
}
int findlast(vector<int>& nums, int start, int end, int target) {
	if (start > end)return -1;
	int mid = (start + end) / 2;
	//if (start == end && nums[mid] == target)return start; else return -1;
	if (nums[mid] > target)return findlast(nums, start, mid - 1, target);
	else if (nums[mid]< target)return findlast(nums, mid + 1, end, target);
	else if (mid + 1<nums.size() && nums[mid + 1] == target)return findlast(nums, mid + 1, end, target);
	else return mid;
}
vector<int> searchRange(vector<int>& nums, int target) {
	vector<int> r;
	r.push_back(findfirst(nums, 0, nums.size() - 1, target));
	r.push_back(findlast(nums, 0, nums.size() - 1, target));
	return r;
}
int SearchInRotatedArray(vector<int>& nums, int target) {
	int lo = 0, hi = nums.size() - 1;
	while (lo < hi) {
		int mid = (lo + hi) / 2;
		if (nums[mid] > nums[hi])lo = mid + 1;
		else hi = mid;
	}
	int rot = lo;
	lo = 0, hi = nums.size() - 1;
	while (lo <= hi) {
		int mid = (lo + hi) / 2;
		int realmid = (mid + rot) % nums.size();
		if (nums[realmid] < target)lo = mid + 1;
		else if (nums[realmid] > target)hi = mid - 1;
		else return realmid;
	}
	return -1;
}
int searchInsert(vector<int>& nums, int target) {
	int lo = 0, hi = nums.size() - 1;
	while (lo <= hi) {
		int mid = (lo + hi) / 2;
		if (nums[mid] < target)lo = mid + 1;
		else if (nums[mid] > target)hi = mid - 1;
		else return mid;
	}
	return lo;
}
bool isValidSudoku(vector<vector<char>>& board) {
	int validR[9], validC[9];
	for (int i = 0; i < 9; i++) {
		memset(validR, 0, sizeof(int) * 9);
		memset(validC, 0, sizeof(int) * 9);
		for (int j = 0; j < 9; j++) {
			if (board[i][j] != '.') {
				if (validR[board[i][j] - '1'])
					return false;
				else validR[board[i][j] - '1'] = 1;
			}
			if (board[j][i] != '.') {
				if (validC[board[j][i] - '1'])
					return false;
				else validC[board[j][i] - '1'] = 1;
			}
		}
	}
	for (int i = 0; i<3; i++)
		for (int j = 0; j < 3; j++) {
			memset(validC, 0, sizeof(int) * 9);
			for (int k = 0; k < 9; k++) {
				if (board[3 * i + k / 3][3 * j + k % 3] != '.') {
					if (validC[board[3 * i + k / 3][3 * j + k % 3] - '1'])
						return false;
					else validC[board[3 * i + k / 3][3 * j + k % 3] - '1'] = 1;
				}
			}
		}
	return 1;
}
vector<vector<int>> permute(vector<int>& nums) {
	vector<vector<int>> r;
	if (nums.empty())return r;
	sort(nums.begin(), nums.end());
	do {
		r.push_back(nums);
	} while (next_permutation(nums.begin(), nums.end()));
	return r;
}
vector<vector<int>> permuteUnique(vector<int>& nums) {
	vector<vector<int>> r;
	if (nums.empty())return r;
	set<vector<int>> set1;
	sort(nums.begin(), nums.end());
	do {
		set1.insert(nums);
	} while (next_permutation(nums.begin(), nums.end()));
	for (auto i = set1.begin(); i != set1.end(); i++)
		r.push_back(*i);
	return r;
}
int firstMissingPositive(vector<int>& nums) {
	int n = nums.size();
	for (int i = 0; i<n; i++)
		if (nums[i] <= 0)
			nums[i] = INT_MAX - 1;
	for (int i = 0; i<n; i++) {
		if (abs(nums[i]) <= n && /*防止重复*/nums[abs(nums[i]) - 1] >0)
			nums[abs(nums[i]) - 1] = -nums[abs(nums[i]) - 1];
	}
	for (int i = 0; i<n; i++)
		if (nums[i]>0)
			return i + 1;
	return n + 1;
}
double myPow(double x, int n) {
	if (n == 0)return 1.0;
	if (n<0)return 1 / x*myPow(1 / x, -n - 1);//to ensure -n<INT_MAX
	return n % 2 == 0 ? myPow(x*x, n / 2) : x*myPow(x*x, n / 2);
}
int trap(vector<int>& height)
{
	int ans = 0, current = 0;
	stack<int> st;
	while (current < height.size()) {
		while (!st.empty() && height[current] > height[st.top()]) {
			int top = st.top();
			st.pop();
			if (st.empty())
				break;
			int distance = current - st.top() - 1;
			int bounded_height = min(height[current], height[st.top()]) - height[top];
			ans += distance * bounded_height;
		}
		st.push(current++);
	}
	return ans;
}
void rotate(vector<vector<int>>& matrix) {
	int n = matrix.size();
	if (n <= 1)return;
	for (int i = 0; i<n / 2; i++)
		for (int j = 0; j < n - 2 * i - 1; j++) {
			int t = matrix[i + j][n - 1 - i];
			matrix[i + j][n - 1 - i] = matrix[i][i + j];
			int t2 = matrix[n - 1 - i][n - 1 - i - j];
			matrix[n - 1 - i][n - 1 - i - j] = t;
			t = matrix[n - 1 - i - j][i];
			matrix[n - 1 - i - j][i] = t2;
			matrix[i][i + j] = t;
		}
	return;
}
string bigint_add(string num1, string num2) {
	string ret;
	int m;
	int up = 0;
	while (num1.length() != 0 && num2.length() != 0) {
		m = (num1[num1.length() - 1] - '0') + (num2[num2.length() - 1] - '0') + up;
		up = m / 10;
		m = m % 10;
		num1.erase(num1.length() - 1); num2.erase(num2.length() - 1);
		ret = to_string(m) + ret;
	}
	while (num1.length() != 0) {
		m = (num1[num1.length() - 1] - '0') + up;
		up = m / 10;
		m = m % 10;
		num1.erase(num1.length() - 1);
		ret = to_string(m) + ret;
	}
	while (num2.length() != 0) {
		m = (num2[num2.length() - 1] - '0') + up;
		up = m / 10;
		m = m % 10;
		num2.erase(num2.length() - 1);
		ret = to_string(m) + ret;
	}
	if (up != 0)ret = to_string(up) + ret;
	while (ret.size() > 1 && ret[0] == '0')
		ret.erase(ret.begin());
	return ret;
}
string bigint_multiply(string num1, string num2) {
	string ret = "0";
	int m;
	int up = 0;
	for (int i = num2.size() - 1; i >= 0; i--) {
		string tmp = "";
		int m;
		int mul = num2[i] - '0';
		int up = 0;
		for (int j = num1.size() - 1; j >= 0; j--) {
			m = (num1[j] - '0')*mul + up;
			up = m / 10;
			m = m % 10;
			tmp = to_string(m) + tmp;
		}
		if (up != 0)tmp = to_string(up) + tmp;
		for (int j = 1; j <= num2.size() - 1 - i; j++)tmp = tmp + "0";
		ret = bigint_add(tmp, ret);
	}
	while (ret.size() > 1 && ret[0] == '0')
		ret.erase(ret.begin());
	return ret;
}
bool isMatch2(string s, string p) {
	if (s.empty() && p.empty())return true;
	if (!s.empty() && p.empty())return false;
	if (s.empty())
		if (p[0] == '*')return isMatch2(s, p.substr(1, p.length() - 1));
		else return false;
	else
		if (p[0] == s[0] || p[0] == '?') {
			p.erase(p.begin()); s.erase(s.begin());
			return isMatch2(s, p);
		}
		else if (p[0] == '*') {
			//while (p.size()>1 && p[1] == '*')p.erase(p.begin());
			return isMatch2(s, p.substr(1, p.length() - 1)) || isMatch2(s.substr(1, s.length() - 1), p);
		}
		else
			return false;
}
bool isMatchttt(string s, string p) {
	int i = 0;
	while (i < p.size()) {
		while (p.size() > 1 && p[i] == '*' &&  p[i + 1] == '*')p.erase(p.begin() + i + 1);
		i++;
	}
	return isMatch2(s, p);
}
int jump(vector<int>& nums) {
	int n = nums.size();
	if (n<2)return 0;
	int level = 0, currentMax = 0, i = 0, nextMax = 0;

	while (currentMax - i + 1>0) {		//nodes count of current level>0
		level++;
		for (; i <= currentMax; i++) {	//traverse current level , and update the max reach of next level
			nextMax = max(nextMax, nums[i] + i);
			if (nextMax >= n - 1)return level;   // if last element is in level+1,  then the min jump=level 
		}
		currentMax = nextMax;
	}
	return 0;
}
int jump2(vector<int>& nums) {
	vector<int> dp(nums.size(), INT_MAX);
	dp[0] = 0;
	for (int i = 1; i < nums.size(); i++)
		for (int j = 0; j <= i - 1; j++) {
			if (nums[j] + j >= i)dp[i] = min(dp[i], dp[j] + 1);
		}
	return dp[nums.size() - 1];
}
int LeetCode62_uniquePaths(int m, int n) {
	//动规
	vector<vector<int>> dp(m, vector<int>(n, 0));
	for (int i = 0; i < m; i++) dp[i][0] = 1;
	for (int i = 0; i < n; i++) dp[0][i] = 1;
	for (int i = 1; i < m; i++)
		for (int j = 1; j < n; j++)
			dp[i][j] = dp[i - 1][j] + dp[i][j - 1];
	return dp[m - 1][n - 1];
}
int LeetCode63_uniquePathsWithObstacles(vector<vector<int>>& obstacleGrid) {
	//动规,有障碍物的地方，不能从那个方向过来
	int m = obstacleGrid.size(), n = obstacleGrid[0].size();
	vector<vector<int>> dp(m, vector<int>(n, 0));
	dp[0][0] = 1;
	for (int i = 1; i < m; i++)
		if (!obstacleGrid[i - 1][0])dp[i][0] = dp[i - 1][0];
	for (int i = 1; i < n; i++)
		if (!obstacleGrid[0][i - 1])dp[0][i] = dp[0][i - 1];
	for (int i = 1; i < m; i++)
		for (int j = 1; j < n; j++) {
			if (!obstacleGrid[i - 1][j])
				dp[i][j] += dp[i - 1][j];
			if (!obstacleGrid[i][j - 1])
				dp[i][j] += dp[i][j - 1];
		}
	if (!obstacleGrid[m - 1][n - 1])
		return dp[m - 1][n - 1];
	else
		//终点被堵住的情况
		return 0;
}
int LeetCode64_0minPathSum(vector<vector<int>>& grid) {
	//动规
	int m = grid.size(), n = grid[0].size();
	vector<vector<int>> dp(m, vector<int>(n, 0));
	dp[0][0] = grid[0][0];
	for (int i = 1; i < m; i++)
		dp[i][0] = dp[i - 1][0] + grid[i][0];
	for (int i = 1; i < n; i++)
		dp[0][i] = dp[0][i - 1] + grid[0][i];
	for (int i = 1; i < m; i++)
		for (int j = 1; j < n; j++)
			dp[i][j] = min(dp[i - 1][j], dp[i][j - 1]) + grid[i][j];
	return dp[m - 1][n - 1];
}
bool LeetCode65_isNumber(string s) {
	char str[256];
	while (!s.empty() && s[0] == ' ')s.erase(s.begin());
	while (!s.empty() && s[s.length() - 1] == ' ')s.erase(s.end() - 1);
	strcpy(str, s.c_str());
	return isNumeric2(str);
}
vector<int> LeetCode66_plusOne(vector<int>& digits) {
	int up = 0;
	int i = digits.size() - 1;
	while (1) {
		if (i < 0)break;
		digits[i] += +up;
		if (i == digits.size() - 1)digits[i]++;
		up = 0;
		if (digits[i] >= 10) {
			up = digits[i] / 10;
			digits[i] %= 10;
		}
		i--;
		if (!up)return digits;
	}
	digits.insert(digits.begin(), up);
	return digits;
}
int LeetCode58_lengthOfLastWord(string s) {
	while (!s.empty() && s[s.length() - 1] == ' ')s.erase(s.end() - 1);
	if (s.empty())return 0;
	int i = s.length() - 1;
	while (i >= 0 && s[i] != ' ')i--;
	return s.length() - 1 - i;
}
int LeetCode28_strStr(string haystack, string needle) {
	if (needle.empty())return 0;
	int i = 0;
	while (i < haystack.length()) {
		//用i从左到右扫描haystack，先找到其中needle的第一个字符
		while (i < haystack.length() && haystack[i] != needle[0])i++;
		//余下的部分不够needle长，肯定不行
		if (i + needle.length() > haystack.length())return -1;
		//逐字符比对
		int j = i;
		for (; j < i + needle.length(); j++)
			if (haystack[j] != needle[j - i])
				break;
		//比对成功
		if (j == i + needle.length())return i;
		//上一轮比对不成功，i+1继续
		i++;
	}
	return -1;
}
int LeetCode70_climbStairs(int n) {
	vector<int>dp(n + 1, 0);
	dp[1] = 1; dp[2] = 2;
	for (int i = 3; i <= n; i++)dp[i] = dp[i - 1] + dp[i - 2];
	return dp[n];
}
void dfs_findlcs(vector<vector<string>>& direction, int i, int j,
	vector<pair<int, int>> curlcs, vector<vector<pair<int, int>>>& lcss) {
	string bp = direction[i][j];
	if (bp == "xd")
		int c = 1;
	if (i == 0 || j == 0) {
		lcss.push_back(curlcs);
		//curlcs.erase(curlcs.end());
		return;
	}
	for (int k = 0; k < direction[i][j].size(); k++) {
		switch (direction[i][j][k]) {
		case 'd':
			dfs_findlcs(direction, i - 1, j, curlcs, lcss);
			break;
		case 'x':
			curlcs.push_back(pair<int, int>(i, j));
			dfs_findlcs(direction, i - 1, j - 1, curlcs, lcss);
			curlcs.pop_back();
			break;
		case 'r':
			dfs_findlcs(direction, i, j - 1, curlcs, lcss);
			break;
		}
	}
	//curlcs.erase(curlcs.end());
}
int LeetCode72_minDistance(string word1, string word2) {
	//1.用dp先找出所有最长公共子序列LCS
	//2.对每个LCS求最少变化步数
	//3.找出最少的变化

	//dp过程
	vector<vector<int>> dp(word1.size() + 1, vector<int>(word2.size() + 1, 0));
	vector<vector<string>> direction(word1.size() + 1, vector<string>(word2.size() + 1, ""));
	for (int i = 1; i <= word1.size(); i++)
		for (int j = 1; j <= word2.size(); j++) {
			if (word1[i - 1] == word2[j - 1]) {//两个新字母相同
				dp[i][j] = dp[i - 1][j - 1] + 1;
				direction[i][j] = "x";
				//direction[i-1][j]>direction[i-1][j-1],说明有另外一条子序列路径
				if (dp[i - 1][j] == dp[i][j])
					direction[i][j] += "d";
				if (dp[i][j - 1] == dp[i][j])
					direction[i][j] += "r";
			}
			else {//新字母不同
				if (dp[i][j - 1] > dp[i - 1][j]) {//左侧路径较大
					dp[i][j] = dp[i][j - 1];
					direction[i][j] = "r";
				}
				else if (dp[i][j - 1] < dp[i - 1][j]) {//上侧路径较大
					dp[i][j] = dp[i - 1][j];
					direction[i][j] = "d";
				}
				//dp[i][j - 1] == dp[i - 1][j],分左上角元素是否小于左和上的情况，若小于，说明左和上是两条路径
				else if (dp[i - 1][j - 1]<dp[i][j - 1]) {
					dp[i][j] = dp[i - 1][j];
					direction[i][j] = "rd";
				}
				else {//若左上，上和左路径数字均相同，说明是同一路径，使用右侧路径
					dp[i][j] = dp[i - 1][j];
					direction[i][j] = "r";
				}
			}
		}

	cout << "dp matrix" << endl;
	for (int i = 0; i <= dp.size() - 1; i++) {
		copy(dp[i].begin(), dp[i].end(), ostream_iterator<int>(cout, " "));
		cout << endl;
	}
	cout << endl;
	cout << "direct matrix" << endl;
	for (int i = 0; i <= direction.size() - 1; i++) {
		copy(direction[i].begin(), direction[i].end(), ostream_iterator<string>(cout, " "));
		cout << endl;
	}


	//找出子序列过程,word1[lcss[k][l].first]==word2[lcss[k][l].second]
	vector<vector<pair<int, int>>> lcss;
	vector<pair<int, int>> curlcs;
	dfs_findlcs(direction, word1.size(), word2.size(), curlcs, lcss);

	cout << endl;
	for (int i = 0; i < lcss.size(); i++) {
		for (int j = lcss[i].size()-1; j >=0 ; j--)
			cout << word1[lcss[i][j].first-1] << " " << lcss[i][j].first-1 << " " << lcss[i][j].second-1 << endl;
		cout << endl;
	}
	
	//对每个子序列进行变换，找出变换最少的步数
	int mincount = INT_MAX;//统计最少步数
	for (int i = 0; i < lcss.size(); i++) {
		int count = 0;//统计步数
		//之前算出来的序号都是加1的，这里做减一处理
		for (int j = 0; j <= lcss[i].size() - 1; j++) {
			lcss[i][j].first--;
			lcss[i][j].second--;
		}
		//1，在两字符串的头尾分别添加虚拟字符，并作为为新的一对
		lcss[i].insert(lcss[i].end(), pair<int, int>(-1, -1));
		lcss[i].insert(lcss[i].begin(), pair<int, int>(word1.length(), word2.length()));
		//don't fucking know why the lcss[i] is reversed, so just fucking reverse it...
		reverse(lcss[i].begin(), lcss[i].end());
		for (int j = 0; j <= lcss[i].size() - 2; j++) {
			//两个字符串对应区间分别的长度
			int l1 = lcss[i][j + 1].first - lcss[i][j].first - 1;
			int l2 = lcss[i][j + 1].second - lcss[i][j].second-1;
			//若l1长，则应删l1-l2,改l2
			//若l2>=l1，则应增l2-l1，改l1
			count += max(l1, l2);
		}
		if (count < mincount)mincount = count;
	}

	return mincount;
}
int LeetCode72_minDistance2(string word1, string word2) {
	//动规实现
	vector<vector<int>> dp(word1.size() + 1, vector<int>(word2.size() + 1, 0));
	//字符串变空串需要删掉其长度个字符
	for (int i = 0; i <= word1.size(); i++)
		dp[i][0] = i;
	for (int i = 0; i <= word2.size(); i++)
		dp[0][i] = i;
	for (int i = 1; i <= word1.size(); i++)
		for (int j = 1; j <= word2.size(); j++) {
			if (word1[i-1] == word2[j-1])//注意下标，第i次循环对应word1[i-1]
				//若最后一对字符相等，则为以下的最小值
				//1，i-1串转为j-1串
				//2，i-1串先转为j串，再把i删了
				//3，i串先转为j-1串，再添加一个j字符
				dp[i][j] = 1 + min(min(dp[i - 1][j - 1] - 1, dp[i - 1][j]), dp[i][j - 1]);
			else
				//否则上述的1，还要加上把i字符转为j字符的1
				dp[i][j] =1+ min(min(dp[i - 1][j - 1], dp[i - 1][j]), dp[i][j - 1]);
		}
	return dp[word1.size()][word2.size()];
}
string LeetCode67_addBinary(string a, string b) {
	string ret;
	int i = a.size() - 1, j = b.size() - 1;
	int up = 0;
	while (i >= 0 && j >= 0) {
		int t = a[i] - '0' + b[j] - '0' + up;
		up = t / 2;
		t = t % 2;
		ret.insert(ret.begin(), t + '0');
		i--;
		j--;
	}
	while (i >= 0) {
		int t = a[i] - '0' +up;
		up = t / 2;
		t = t % 2;
		ret.insert(ret.begin(), t + '0');
		i--;
	}
	while (j >= 0) {
		int t = b[j] - '0' + up;
		up = t / 2;
		t = t % 2;
		ret.insert(ret.begin(), t + '0');
		j--;
	}
	if (up > 0) {
		ret.insert(ret.begin(), up + '0');
	}
	return ret;
}
int LeetCode69_mySqrt(int x) {
	if (x == 1) return 1;
	int lo = 1, hi = x / 2 + 1;
	while (lo < hi)
	{
		int mi = (hi + lo) >> 1;
		x / mi < mi ? 
			hi = mi : 
			lo = mi + 1;
	}
	return --lo;
}
string LeetCode71_simplifyPath(string path) {
	list<string> queue;
	string ret;
	int i = 1,j=1;
	while ((j=path.find('/', i))!=string::npos) {
		queue.push_back(path.substr(i, j-i));
		i = j + 1;
	}
	queue.push_back(path.substr(i, path.length() - i));//最后一个项之后可能没/

	list<string>::iterator t;
	for (auto i = queue.begin(); i != queue.end();/*操作内进行迭代器步进*/) {
		if (*i == ".." ){
			if (i != queue.begin()) {
				t = i;
				queue.erase(--i);
				i = t;
				queue.erase(i++);//删除会使迭代器失效，这里也可用i=queue.erase(i);
				//prev = i;
				continue;
			}
			else {
				queue.erase(i++);
				//prev = i;
				continue;
			}
		}
		if (*i == "." || (*i).empty()) {
			queue.erase(i++);
			//prev = i;
			continue;
		}
		//prev = i;
		++i;
	}
	if (!queue.empty()) {
		for (string s : queue)
			ret += "/" + s;
		return ret;
	}
	else
		return "/";
}
void LeetCode73_setZeroes(vector<vector<int>>& matrix) {
	int m = matrix.size(), n = matrix[0].size();
	vector<bool>row(m, 0);
	vector<bool>col(n, 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n;j++)
			if (matrix[i][j]==0) {
				row[i] = true;
				col[j] = true;
			}
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) 
			if (row[i] || col[j])matrix[i][j] = 0;
}
bool LeetCode74_searchMatrix(vector<vector<int>>& matrix, int target) {
	int m = matrix.size();
	if (m == 0)return false;
	int n = matrix[0].size();
	int i = 0,j=m*n-1;
	int mid;
	while (i <= j) {
		mid = (i + j) >> 1;
		if (matrix[mid / n][mid%n] == target)
			return true;
		else if (matrix[mid / n][mid%n] < target)
			i= mid + 1;
		else
			j= mid - 1;
	}
	return false;
}
void LeetCode75_sortColors(vector<int>& nums) {
	int i = 0, j = nums.size()-1;
	//int ones = 0;
	while (j>i) {
		while (j>i && nums[i] == 0){
			/*if(nums[i]==1){
				nums[i] = 0;
				ones++;}*/
			i++;}
		while (j>i && nums[j] != 0){
			/*if (nums[j] == 1) {
				nums[j] = 2;
				ones++;
			}*/
			j--;}
		swap(nums[i], nums[j]);
	}
	if (nums[i] == 0)i = i + 1;
	j = nums.size() - 1;
	while (j>i) {
		while (j>i && nums[i] == 1) {
			/*if(nums[i]==1){
			nums[i] = 0;
			ones++;}*/
			i++;
		}
		while (j>i && nums[j] != 1) {
			/*if (nums[j] == 1) {
			nums[j] = 2;
			ones++;
			}*/
			j--;
		}
		swap(nums[i], nums[j]);
	}
	//for (int k = i+1;)
	int c = 1;
	return;
}
void LeetCode77_dfs(int cur, int lv,int n, vector<int>& arr, vector<vector<int>>& ret) {
	if (lv == 1){
		ret.push_back(arr);
		//arr.pop_back();
		return;}
	for (int i = cur+1; i <= n; i++) {
		arr.push_back(i);
		LeetCode77_dfs(i, lv - 1, n, arr, ret);
		arr.pop_back();
	}
}
vector<vector<int>> LeetCode77_combine(int n, int k) {
	vector<vector<int>>ret;
	vector<int> arr;
	LeetCode77_dfs(0, k + 1, n, arr, ret);
	return ret;
}
void LeetCode78_dfs(int cur, int lv, vector<int>& arr, vector<int>& nums, vector<vector<int>>& ret) {
	if (lv == 1) {
		ret.push_back(arr);
		//arr.pop_back();
		return;
	}
	for (int i = cur + 1; i < nums.size(); i++) {
		arr.push_back(nums[i]);
		LeetCode78_dfs(i, lv - 1,  arr,nums, ret);
		arr.pop_back();
	}
}
vector<vector<int>> LeetCode78_subsets(vector<int>& nums) {
	vector<vector<int>> ret;
	vector<int> arr;
	ret.push_back(arr);
	for (int i = 1; i <= nums.size(); i++)
		LeetCode78_dfs(-1, i+1, arr, nums, ret);
	return ret;
}
int LeetCode80_removeDuplicates(vector<int>& nums) {
	int i = 0;
	for (int n : nums)
		if (i < 2 || n > nums[i - 2])
			nums[i++] = n;
	return i;
}
bool LeetCode81_search(vector<int>& nums, int target) {
	if (!nums.size()) return false;
	int left = 0;
	int right = nums.size() - 1;
	if (nums[left] == target || nums[right] == target) return true;
	while (left <= right) {
		int mid = (left + right) >> 1;

		if (nums[mid] == target) return true;

		if (nums[left] <= target && nums[mid] > target)  // this would be a normal situation...
			right = mid - 1;
		else if (nums[right] >= target && nums[mid] < target) 
			left = mid + 1;
		else if (nums[right] >= target && nums[mid] > target) 
			left++;
		else
			right--;
	}
	return false;
}
ListNode* LeetCode83_deleteDuplicates(ListNode* head) {
	ListNode* fhead = new ListNode(INT_MIN);
	fhead->next = head;
	ListNode* p=head, *prev=fhead;
	while (p) {
		if (p->val == prev->val) 
			prev->next = p->next;
		else
			prev = p;
		p = p->next;
	}
	return fhead->next;
}
ListNode* LeetCode82_deleteDuplicates(ListNode* head) {
	ListNode* fhead = new ListNode(INT_MIN);
	fhead->next = head;
	ListNode* p = head, *prev = fhead;
	while (p) {
		if (p->next && p->val == p->next->val){
			while (p->next && p->val == p->next->val)
				p->next = p->next->next;
			prev->next = p->next;}
		else
			prev = p;
		p = p->next;
	}
	return fhead->next;
}
int LeetCode84_largestRectangleArea(vector<int>& heights) {
	int maxArea = 0;
	stack<int> stack;

	int i = 0;
	while (i < heights.size() || !stack.empty()) {
		//一路将递增序列推入栈。
		if (stack.empty() || (i < heights.size() && heights[i] >= heights[stack.top()])) {
			stack.push(i++);
		}
		else {
			//出栈
			int top = stack.top();
			stack.pop();
			//出栈后，栈顶元素实际是上一个比它矮的元素
			//这时，[stack.top()+1,i-1]是以heights[top]为高能延展的最大区间
			int width = (stack.empty()) ? i : i - 1 - stack.top();
			int area = width * heights[top];

			maxArea =max(area, maxArea);
		}
	}
	return maxArea;
}
void LeetCode88_merge(vector<int>& nums1, int m, vector<int>& nums2, int n) {
	int index = m + n - 1,i=m-1,j=n-1;
	while (index >= 0) {
		if (i >= 0 && j>=0) {
			if (nums1[i] > nums2[j])
				nums1[index--] = nums1[i--];
			else
				nums1[index--] = nums2[j--];
		}
		else if(i>=0)
			nums1[index--] = nums1[i--];
		else
			nums1[index--] = nums2[j--];
	}
}
void OutPutList(ListNode* h) {
	while (h->next) {
		cout << h->val << "->";
		h = h->next;
	}
	cout << h->val << endl;
}
ListNode* LeetCode86_partition(ListNode* head, int x) {
	ListNode* big = new ListNode(INT_MIN), *small = new ListNode(INT_MIN),*p=head,
		*bp=big,*sp=small;
	while (p) {
		if (p->val < x)
			bp->next = p,
			p = p->next,
			bp = bp->next;
		else
			sp->next = p,
			p = p->next,
			sp = sp->next;
	}
	bp->next = NULL, sp->next = NULL;
	//OutPutList(big), OutPutList(small);
	/*bp = big->next;
	while (bp->next)
		bp = bp->next;*/
	bp->next = small->next;
	return big->next;
}
ListNode* LeetCode92_reverseBetween(ListNode* head, int m, int n) {
	if (m == n)return head;
	ListNode* m_prev=NULL, *n_next;
	ListNode* p, *pp, *ppp;
	//如果不是从链表头开始处理
	if (m != 1) {
		m_prev = head;
		for (int i = 1; i <= m - 2; i++)
			m_prev = m_prev->next;
		p = m_prev->next;
	}
	else
		p = head;
	pp = p->next;
	ppp = pp->next;
	//m到n一共n-m+1个数，需做n-m次
	for (int i = 1; i <=n - m ; i++) {
		pp->next = p;
		p = pp;
		pp = ppp;
		if (ppp)
			ppp = ppp->next;
	}
	n_next = pp;
	//如果不是从链表头开始处理
	if (m_prev) {
		pp = m_prev->next;
		m_prev->next = p;//p此时指向反转的最后一个
	}
	//否则
	else{
		pp = head;
		head = p;}
	pp->next = n_next;//用pp储存反转的第一个
	return head;
}
vector<string> LeetCode93_restoreIpAddresses(string s) {
	vector<string> ret;
	int dot1 = 1, dot2 = 2, dot3 = 3;
	//sb测试用例
	if (s.length()<4 || s.length()>12)
		return ret;
	while(1){
		//不能是数字序列开头为0（且不是单独一个0）
		if(!(dot1>0+1 && s[0]=='0' ||
			dot2>dot1 + 1 && s[dot1] == '0' || 
			dot3>dot2 + 1 && s[dot2] == '0' || 
			s.length()>dot3 + 1 && s[dot3] == '0')){
			//在dot位置插入点，并转换为整数
			int num1 = atoi(s.substr(0, dot1).c_str());
			int num2 = atoi(s.substr(dot1, dot2-dot1).c_str());
			int num3 = atoi(s.substr(dot2, dot3 - dot2).c_str());
			int num4 = atoi(s.substr(dot3, s.length()- dot3).c_str());
			if (num1 >= 0 && num1 <= 255 &&
				num2 >= 0 && num2 <= 255 &&
				num3 >= 0 && num3 <= 255 &&
				num4 >= 0 && num4 <= 255 ){
				string t = s;
				t.insert(dot3,".");
				t.insert(dot2, ".");
				t.insert(dot1, ".");
				ret.push_back(t);
			}}
		//生成下一个dot序列
		if (dot3 < s.length() - 1)
			dot3++;
		else if (dot2 < s.length() - 2) {
			dot2++;
			dot3 = dot2 + 1;
		}
		else if (dot1 < s.length() - 3) {
			dot1++;
			dot2 = dot1 + 1;
			dot3 = dot2 + 1;
		}
		else
			break;
		int c = 1;
	}
	return ret;
}
void LeetCode94_inorderTraversal(TreeNode* p,vector<int>& ret) {
	if (!p)
		return;
	if (p->left)
		LeetCode94_inorderTraversal(p->left,ret);
	ret.push_back(p->val);
	if (p->right)
		LeetCode94_inorderTraversal(p->right,ret);
}
vector<int> LeetCode94_inorderTraversal(TreeNode* root) {
	vector<int> ret;
	LeetCode94_inorderTraversal(root, ret);
	return ret;
}
int LeetCode96_numTrees(int n) {
	vector<int> BSTs(n+1, 0);
	BSTs[0] = 1;
	BSTs[1] = 1;
	for (int i = 2; i <= n; i++) {
		int count = 0;
		for (int j = 1; j <= i; j++)
			count += BSTs[j - 1] * BSTs[i-j];
		BSTs[i] = count;
	}
	return BSTs[n];
}
TreeNode* LeetCode95_CloneAndChangeNumber(TreeNode* p, int start) {
	if (!p)
		return NULL;
	TreeNode* t = new TreeNode(p->val-1+start);
	t->left = LeetCode95_CloneAndChangeNumber(p->left, start);
	t->right= LeetCode95_CloneAndChangeNumber(p->right, start);
	return t;
}
vector<TreeNode*> LeetCode95_generateTrees(int n) {
	vector<vector<TreeNode*>> BSTs;
	vector<TreeNode*> tmp;
	//构造0和1的情况
	tmp.push_back(NULL);
	BSTs.push_back(tmp);
	tmp.clear();
	tmp.push_back(new TreeNode(1));
	BSTs.push_back(tmp);
	//逐一构造2到n的情况
	for (int i = 2; i <= n; i++) {
		tmp.clear();
		//对应每个i，要将BSTs[j-1]和BSTs[i-j]的树逐一相乘
		for (int j = 1; j <= i; j++) {
			for (int k = 0; k < BSTs[j-1].size(); k++) {
				for (int l = 0; l < BSTs[i-j].size(); l++){
					TreeNode* h = new TreeNode(j);
					h->left = LeetCode95_CloneAndChangeNumber(BSTs[j - 1][k],1);
					h->right= LeetCode95_CloneAndChangeNumber(BSTs[i - j][l], j+1);
					tmp.push_back(h);
				}}}	
		BSTs.push_back(tmp);
	}
	BSTs[0].clear();
	return BSTs[n];
}
bool LeetCode100_isSameTree(TreeNode* p, TreeNode* q) {
	if (p && q && p->val == q->val)
		return LeetCode100_isSameTree(p->left, q->left) && LeetCode100_isSameTree(p->right, q->right);
	else if (!p && !q)
		return true;
	else
		return false;
}
int LeetCode98_minofBST(TreeNode* root) {
	while (root->left)
		root = root->left;
	return root->val;
}
int LeetCode98_maxofBST(TreeNode* root) {
	while (root->right)
		root = root->right;
	return root->val;
}
bool LeetCode98_isValidBST(TreeNode* root) {
	if (!root)
		return true;
	if (root->left)
		if (root->left->val >= root->val || !LeetCode98_isValidBST(root->left) || LeetCode98_maxofBST(root->left) >= root->val)
			return false;
	if (root->right)
		if (root->right->val <= root->val || !LeetCode98_isValidBST(root->right) || LeetCode98_minofBST(root->right) <= root->val)
			return false;
	return true;
}
bool LeetCode97_isInterleave(string s1, string s2, string s3) {
	if (s3.length() != s1.length() + s2.length()) {
		return false;
	}
	vector<vector<bool>> dp(s1.length() + 1,vector<bool>(s2.length() + 1,false));
	for (int i = 0; i <= s1.length(); i++) {
		for (int j = 0; j <= s2.length(); j++) {
			if (i == 0 && j == 0) {
				dp[i][j] = true;
			}
			else if (i == 0) {
				dp[i][j] = dp[i][j - 1] && s2[j - 1]== s3[i + j - 1];
			}
			else if (j == 0) {
				dp[i][j] = dp[i - 1][j] && s1[i - 1 ]== s3[i + j - 1];
			}
			else {
				dp[i][j] = (dp[i - 1][j] && s1[i - 1] == s3[i + j - 1]) || (dp[i][j - 1] && s2[j - 1] == s3[i + j - 1]);
			}
		}
	}
	return dp[s1.length()][s2.length()];
}
TreeNode* LeetCode101_flip(TreeNode* root) {
	if (!root)return nullptr;
	TreeNode* r = new TreeNode(root->val);
	if (root->right)
		r->left = LeetCode101_flip(root->right);
	if (root->left)
		r->right= LeetCode101_flip(root->left);
	return r;
}
bool LeetCode101_isSymmetric(TreeNode* root) {
	TreeNode* r = LeetCode101_flip(root);
	return LeetCode100_isSameTree(r, root);
}
void LeetCode102_dfs(TreeNode*p,int lv, vector<vector<int>>& ret) {
	if (!p)return;
	if (lv >= ret.size())
		ret.resize(ret.size() + 1);
	ret[lv].push_back(p->val);
	LeetCode102_dfs(p->left, lv + 1, ret);
	LeetCode102_dfs(p->right, lv + 1, ret);
}
vector<vector<int>> LeetCode102_levelOrder(TreeNode* root) {
	vector<vector<int>>ret;
	LeetCode102_dfs(root, 0, ret);
	return ret;
}
vector<vector<int>> LeetCode103_zigzagLevelOrder(TreeNode* root) {
	vector<vector<int>>ret;
	LeetCode102_dfs(root, 0, ret);
	for (int i = 1; i < ret.size(); i += 2)
		reverse(ret[i].begin(), ret[i].end());
	return ret;
}
int LeetCode104_maxDepth(TreeNode* root) {
	if (!root)return 0;
	return max(LeetCode104_maxDepth(root->left), LeetCode104_maxDepth(root->right)) + 1;
}
TreeNode* LeetCode105_buildTree(const vector<int>& preorder, const vector<int>& inorder) {
	if (preorder.empty() || inorder.empty())
		return nullptr;
	TreeNode* h = new TreeNode(preorder[0]);
	int h_index = find(inorder.begin(), inorder.end(), preorder[0])-inorder.begin();
	if(h_index + 1<=preorder.size())
		h->left = LeetCode105_buildTree(vector<int>(preorder.begin() + 1, preorder.begin() +h_index+1),
			vector<int>(inorder.begin(), inorder.begin()+h_index)); 
	if (h_index+1<preorder.size())
		h->right = LeetCode105_buildTree(vector<int>(preorder.begin() +h_index + 1, preorder.end()),
				vector<int>(inorder.begin() + h_index+1, inorder.end()));
	return h;
}
TreeNode* LeetCode106_buildTree(const vector<int>& inorder, const vector<int>& postorder) {
	if (postorder.empty() || inorder.empty())
		return nullptr;
	TreeNode* h = new TreeNode(postorder.back());
	int h_index = find(inorder.begin(), inorder.end(), postorder.back()) - inorder.begin();
	if (h_index <= postorder.size())
		h->left = LeetCode106_buildTree(vector<int>(inorder.begin(), inorder.begin() + h_index),
			vector<int>(postorder.begin(), postorder.begin() + h_index));
	if (h_index + 1<postorder.size())
		h->right = LeetCode106_buildTree(vector<int>(inorder.begin() + h_index + 1, inorder.end()),
			vector<int>(postorder.begin() + h_index, postorder.end() - 1));
	return h;
}
vector<vector<int>> LeetCode107_levelOrderBottom(TreeNode* root) {
	vector<vector<int>>ret;
	LeetCode102_dfs(root, 0, ret);
	reverse(ret.begin(), ret.end());
	return ret;
}
TreeNode* LeetCode108_sortedArrayToBST(const vector<int>& nums) {
	if (nums.empty())
		return nullptr;
	int mid = (nums.size() - 1) / 2;
	TreeNode* h = new TreeNode(nums[mid]);
	h->left = LeetCode108_sortedArrayToBST(vector<int>(nums.begin(), nums.begin() + mid));
	h->right = LeetCode108_sortedArrayToBST(vector<int>(nums.begin()+mid+1, nums.end()));
	return h;
}
TreeNode* LeetCode109_sortedListToBST(ListNode* head) {
	vector<int> nums;
	while (head) {
		nums.push_back(head->val);
		head = head->next;
	}
	return LeetCode108_sortedArrayToBST(nums);
}
bool LeetCode110_isBalanced(TreeNode* root) {
	if (!root)
		return true;
	if (abs(LeetCode104_maxDepth(root->left) - LeetCode104_maxDepth(root->right)) > 1)
		return false;
	return LeetCode110_isBalanced(root->left) && LeetCode110_isBalanced(root->right);
}
int LeetCode111_recursive(TreeNode* root) {
	if (!root)
		return INT_MAX;
	if (!root->right && !root->left)
		return 1;
	return min(LeetCode111_recursive(root->left), LeetCode111_recursive(root->right)) + 1;
}
int LeetCode111_minDepth(TreeNode* root) {
	if (!root)
		return 0;
	else
		return LeetCode111_recursive(root);
}
bool LeetCode112_recursive(TreeNode* root, int sum) {
	//叶子节点
	if (!root -> left && !root->right) {
		if (sum == root->val)
			return true;
		else
			return false;
	}
	else {
		bool t1 = false, t2 = false;
		if (root->left)
			t1 = LeetCode112_recursive(root->left, sum - root->val);
		if (root->right)
			t2 = LeetCode112_recursive(root->right, sum - root->val);
		return t1 || t2;
	}
}
bool LeetCode112_hasPathSum(TreeNode* root, int sum) {
	if (!root)
		return false;
	else
		return LeetCode112_recursive(root, sum);
}
TreeNode* MakeTree(vector<int> arr,int nullval) {
	vector<TreeNode*> p;
	for (int i = 0; i < arr.size(); i++) {
		p.push_back(new TreeNode(arr[i]));
		if (arr[i] != nullval && i!=0) {
			if (i % 2 == 0)
				p[(i - 1) / 2]->right = p[i];
			else
				p[(i - 1) / 2]->left = p[i];
		}
	}
	return p[0];
}
void LeetCode113_recursive(TreeNode* root, int sum, vector<int>&path,vector<vector<int>> &ret) {
	path.push_back(root->val);
	//叶子节点
	if (!root->left && !root->right) {
		if (sum == root->val)
			ret.push_back(path);
	}
	else {
		if (root->left)
			LeetCode113_recursive(root->left, sum - root->val,path,ret);
		if (root->right)
			LeetCode113_recursive(root->right, sum - root->val,path,ret);
	}
	path.pop_back();
}
vector<vector<int>> LeetCode113_pathSum(TreeNode* root, int sum) {
	vector<vector<int>> ret;
	if (!root)
		return ret;
	else{
		vector<int> tmp;
		LeetCode113_recursive(root, sum,tmp,ret);}
	return ret;
}
void LeetCode114_dfs(TreeNode* root,vector<TreeNode*> &q) {
	q.push_back(root);
	if (!root->left && !root->right) {
		return;
	}
	else {
		if (root->left)
			LeetCode114_dfs(root->left,q);
		if (root->right)
			LeetCode114_dfs(root->right,q);
	}
}
void LeetCode114_flatten(TreeNode* root) {
	vector<TreeNode*> q;
	if (!root)return;
	LeetCode114_dfs(root, q);
	int t = q.size() - 2;
	for (int i = 0; i <= t; i++){
		q[i]->left= NULL;
		q[i]->right = q[i + 1];}
	q[q.size() - 1]->left= NULL;
	q[q.size() - 1]->right = NULL;
}
void LeetCode257_recursive(TreeNode* root, vector<int>&path, vector<vector<int>> &ret) {
	path.push_back(root->val);
	//叶子节点
	if (!root->left && !root->right) {
			ret.push_back(path);
	}
	else {
		if (root->left)
			LeetCode257_recursive(root->left, path, ret);
		if (root->right)
			LeetCode257_recursive(root->right, path, ret);
	}
	path.pop_back();
}
vector<string> LeetCode257_binaryTreePaths(TreeNode* root) {
	vector<vector<int>> ret;
	vector<string> r;
	if (!root)
		return r;
	else {
		vector<int> tmp;
		LeetCode257_recursive(root,tmp, ret);
		for (int i = 0; i < ret.size(); i++) {
			string t = to_string(ret[i][0]);
			for (int j = 1; j < ret[i].size(); j++) 
				t += "->" + to_string(ret[i][j]);
			r.push_back(t);
		}
	}
	return r;
}
bool LeetCode141_hasCycle(ListNode *head) {
	if (!head)return false;
	ListNode* f=head, *s=head;
	while (1) {
		f = f->next;
		if (f)
			f = f->next;
		s = s->next;
		if (f&& s) {
			if (s == f)
				return true;
		}
		else
			break;
	}
	return false;
}
ListNode *LeetCode142_detectCycle(ListNode *head) {
	if (!head)return NULL;
	ListNode* f = head, *s = head;
	while (1) {
		f = f->next;
		if (f)
			f = f->next;
		s = s->next;
		if (f&& s) {
			if (s == f)
				break;
		}
		else
			return NULL;
	}
	f = head;
	while (s != f) {
		s = s->next;
		f = f->next;
	}
	return s;
}
vector<int> LeetCode144_preorderTraversal(TreeNode* root) {
	vector<int> ret;
	if (!root)
		return ret;
	stack<TreeNode*>  st;
	st.push(root);
	while (!st.empty()) {
		TreeNode* t = st.top();
		ret.push_back(t->val);
		st.pop();
		if (t->right)
			st.push(t->right);
		if (t->left)
			st.push(t->left);
	}
	return ret;
}
vector<int> LeetCode145_postorderTraversal(TreeNode* root) {
	vector<int> ret;
	if (!root)
		return ret;
	stack<TreeNode*>  st;
	TreeNode* curr=root;
	while (!st.empty() || curr) {
		if (curr) {
			ret.push_back(curr->val);
			st.push(curr->left);
			curr = curr->right;
		}
		else{
			curr = st.top();
			st.pop();}
	}
	reverse(ret.begin(), ret.end());
	return ret;
}
vector<int> preorderTraversal_iteration(TreeNode* root) {
	vector<int> ret;
	if (!root)
		return ret;
	stack<TreeNode*>  st;
	TreeNode* curr = root;
	while (!st.empty() || curr) {
		if (curr) {
			ret.push_back(curr->val);
			st.push(curr->right);
			curr = curr->left;
		}
		else {
			curr = st.top();
			st.pop();
		}
	}
	return ret;
}
vector<int> inorderTraversal_iteration(TreeNode* root) {
	vector<int> ret;
	if (!root)
		return ret;
	stack<TreeNode*>  st;
	TreeNode* curr = root;
	while (!st.empty() || curr) {
		if (curr) {
			st.push(curr);
			curr = curr->left;
		}
		else {
			curr = st.top();
			st.pop();
			ret.push_back(curr->val);
			curr = curr->right;
		}
	}
	return ret;
}
ListNode* LeetCode206_reverseList(ListNode* head) {
	ListNode* p=head, *pp, *ppp=NULL;
	if (!p)
		return NULL;
	pp = p->next;
	if (pp)
		ppp = pp->next;
	while (pp) {
		//ppp = pp->next;
		pp->next = p;
		p = pp;
		pp = ppp;
		if(ppp)
			ppp = ppp->next;
	}
	head->next = NULL;
	return p;
}
void LeetCode143_reorderList(ListNode* head) {
	if (!head)return;
	int n = 0;
	ListNode* p = head;
	while (p) {
		n++;
		p = p->next;
	}
	p = head;
	for (int i = 0; i < n / 2; i++)
		p = p->next;
	ListNode* h2 = p->next;
	p->next = NULL;
	h2=LeetCode206_reverseList(h2);
	if (!h2)return;
	p = head;
	while (1) {
		ListNode* t = p->next;
		p->next = h2;
		p = t;
		t = h2->next;
		h2->next = p;
		if (t)h2 = t;
		else
			break;
	}
	if (p)
		h2->next = p;
}
int LeetCode85_maximalRectangle_1(vector<vector<char> > &matrix) {
	//分行处理，每行用height[i]记载可达高，left[i]表示可向左延伸最远,right[i]表示可向右延伸最远+1
	//最后用max(height[i]*(right[i]-left[i]))得出这一行处理的最大数据
	if (matrix.empty()) return 0;
	const int m = matrix.size();
	const int n = matrix[0].size();
	vector<int> left(n,0), right(n,n), height(n,0);
	int maxA = 0;
	for (int i = 0; i<m; i++) {
		int cur_left = 0, cur_right = n;
		for (int j = 0; j<n; j++) { // compute height (can do this from either side)
			if (matrix[i][j] == '1') height[j]++;
			else height[j] = 0;
		}
		//此时left中存储的是上一行计算的信息
		for (int j = 0; j<n; j++) { // compute left (from left to right)
			if (matrix[i][j] == '1') left[j] = max(left[j], cur_left);
			else { left[j] = 0; cur_left = j + 1; }
		}
		// compute right (from right to left)
		for (int j = n - 1; j >= 0; j--) {
			if (matrix[i][j] == '1') right[j] = min(right[j], cur_right);
			else { right[j] = n; cur_right = j; }
		}
		// compute the area of rectangle (can do this from either side)
		for (int j = 0; j<n; j++)
			maxA = max(maxA, (right[j] - left[j])*height[j]);
	}
	return maxA;
}
int LeetCode85_maximalRectangle_2(vector<vector<char> > &matrix) {
	if (matrix.empty() || matrix[0].empty())
		return 0;
	int n = matrix[0].size();
	vector<int>height(n, 0);
	vector<int> st;
	int ans = 0;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < n; j++)
			height[j] = matrix[i][j] == '1' ? ++height[j] : 0;
		int t = LeetCode84_largestRectangleArea(height);
		ans = max(t, ans);
	}
	return ans;
}
vector<string> LeetCode68_fullJustify(vector<string>& words, int maxWidth) {
	vector<string> ret;
	int i = 0,j=1;
	int len;
	while(i<words.size()){
		len = words[i].size();
		//决定一行，使用词i到j
		while (j<words.size() && len + 1 + words[j].length() <= maxWidth){
			len += words[j].length()+1;
			words[j-1] = words[j-1] + ' ';
			j++;
		}
		//到了最后一行
		if (j >=words.size())
			break;
		//补全空格
		int ii = i;
		for (int k = len; k < maxWidth; k++) {
			words[ii] = words[ii]+ ' ' ;
			++ii;
			if (ii >= j-1)
				ii = i;
		}
		//写出一行
		string t;
		for (int k = i; k <j; k++) {
			t += words[k];
		}
		ret.push_back(t);
		i = j++;
	}
	//处理最后一行，左对齐
	string t;
	for (int k = i; k <j; k++) {
		t += words[k];
	}
	for (int k = len; k < maxWidth; k++) {
		t += ' ';
	}
	ret.push_back(t);
	return ret;
}
vector<char> LeetCode37_ready(int i, int j,vector<vector<char>>& board) {
	vector<bool> b(10,true);
	b[0] = false;
	for (int k = 0; k < 9; k++) {
		if(board[i][k] !='.')
			b[board[i][k] - '0'] = false;
		if (board[k][j] != '.')
			b[board[k][j] - '0'] = false;
	}
	for (int k = 0; k < 3; k++)
		for (int l = 0; l < 3; l++) 
			if(board[i / 3 * 3 + k][j / 3 * 3 + l]!='.')
				b[board[i / 3 * 3 + k][j / 3 * 3 + l] - '0'] = false;
	vector<char> ret;
	for (int k = 1; k <= 9; k++)
		if (b[k])
			ret.push_back(k + '0');
	return ret;
}
bool LeetCode37_track(int n,vector<vector<char>>& board) {
	if (n == 81)
		return true;
	int i = n / 9, j = n % 9;
	if (board[i][j] != '.')
		return LeetCode37_track(n + 1,board);
	vector<char> r = LeetCode37_ready(i, j, board);
	for (char c : r) {
		board[i][j] = c;
		if (LeetCode37_track(n + 1, board))
			return true;
	}
	board[i][j] = '.';
	return false;
}
void LeetCode37_solveSudoku(vector<vector<char>>& board) {
	LeetCode37_track(0, board);
}
bool LeetCode79_recursive(int i,int j,vector<vector<char>>& board, string word, vector<vector<bool>> &mask) {
	if (word.empty())
		return true;
	if (i > 0 && mask[i - 1][j] && board[i - 1][j] == word[0]){
		mask[i - 1][j] = false;
		bool t=LeetCode79_recursive(i - 1, j, board, word.substr(1, word.length() - 1), mask);
		if (t)
			return true;
		else
			mask[i - 1][j] = true;
	}
	if (i < board.size()-1 && mask[i + 1][j] && board[i + 1][j] == word[0]) {
		mask[i+1][j] = false;
		bool t=LeetCode79_recursive(i + 1, j, board, word.substr(1, word.length() - 1), mask);
		if (t)
			return true;
		else
			mask[i + 1][j] = true;
	}
	if (j > 0 && mask[i ][j-1] && board[i ][j-1] == word[0]) {
		mask[i ][j-1] = false;
		bool t= LeetCode79_recursive(i ,j-1, board, word.substr(1, word.length() - 1), mask);
		if (t)
			return true;
		else
			mask[i][j - 1] = true;
	}
	if (j < board[0].size() - 1 && mask[i ][j+1] && board[i ][j+1] == word[0]) {
		mask[i ][j+1] = false;
		bool t= LeetCode79_recursive(i , j + 1, board, word.substr(1, word.length() - 1), mask);
		if (t)
			return true;
		else
			mask[i][j + 1] = true;
	}
	return false;
}
bool LeetCode79_exist(vector<vector<char>>& board, string word) {
	int m = board.size(), n = board[0].size();
	vector<vector<bool>> mask(m, vector<bool>(n, true));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (board[i][j] == word[0]) {
				mask[i][j] = false;
				if(LeetCode79_recursive(i, j , board, word.substr(1, word.length() - 1), mask))
					return true;
				else 
					mask[i][j] = true;
			}}}
	return false;
}
void LeetCode90_dfs(int cur, int lv, vector<int>& arr, vector<int>& nums, set<vector<int>>& ret) {
	if (lv == 1) {
		ret.insert(arr);
		//arr.pop_back();
		return;
	}
	for (int i = cur + 1; i < nums.size(); i++) {
		arr.push_back(nums[i]);
		LeetCode90_dfs(i, lv - 1, arr, nums, ret);
		arr.pop_back();
	}
}
vector<vector<int>> LeetCode90_subsetsWithDup(vector<int>& nums) {
	vector<vector<int>> ret;
	set<vector<int>> set;
	vector<int> arr;
	ret.push_back(arr);
	sort(nums.begin(), nums.end());
	for (int i = 1; i <= nums.size(); i++)
		LeetCode90_dfs(-1, i + 1, arr, nums, set);
	for (auto i = set.begin(); i != set.end(); i++)
		ret.push_back(*i);
	return ret;
}
bool LeetCode76_isdesirable_2(string& t, vector<int> &pos,string& s,int left,int right) {
	vector<bool>flag(t.size(), false);
	for (int i = left; i <= right; ++i) {
		if (binary_search(t.begin(),t.end(),s[pos[i]]))
			for(auto j=lower_bound(t.begin(),t.end(),s[pos[i]]);j!=upper_bound(t.begin(), t.end(), s[pos[i]]);j++)
				if (!flag[j - t.begin()]) {
					flag[j - t.begin()] = true;
					break;
				}
	}
	for (int i = 0; i < flag.size(); i++)
		if (!flag[i])
			return false;
	return true;
}
string LeetCode76_minWindow_2(string s, string t) {
	if (t.empty())
		return"";
	sort(t.begin(), t.end());
	int left = 0, right = 0;
	int minl = -1, minr = -1;
	int min1= INT_MAX;
	vector<int> pos;
	for (int i = 0; i < s.length(); i++)
		if (binary_search(t.begin(),t.end(),s[i]))
			pos.push_back(i);
	if (pos.empty())
		return "";
	while (right <pos.size()) {
		while (right <pos.size() && !LeetCode76_isdesirable_2(t, pos, s, left, right))
			++right;
		if (right >= pos.size())
			break;
		if (min1 > pos[right] - pos[left] + 1) {
			min1 = pos[right] - pos[left] + 1;
			minl = left;
			minr = right;
		}
		while (left+1<=right && LeetCode76_isdesirable_2(t, pos, s, left + 1, right)) 
			++left;
		if (min1 > pos[right] - pos[left] + 1) {
			min1 = pos[right] - pos[left] + 1;
			minl = left;
			minr = right;
		}
		++left;
	}
	if (minl != -1 && minr != -1)
		return s.substr(pos[minl], pos[minr] - pos[minl] + 1);
	else
		return "";
}

bool LeetCode76_isdesirable(multiset<char> set, vector<int> &pos, string& s, int left, int right) {
	for (int i = left; i <= right; ++i) {
		if (set.find(s[pos[i]]) != set.end())
			set.erase(set.lower_bound(s[pos[i]]));
	}
	if (set.empty())
		return true;
	else
		return false;
}
string LeetCode76_minWindow(string s, string t) {
	multiset<char> set;
	for (char c : t)
		set.insert(c);
	int left = 0, right = 0;
	int minl = -1, minr = -1;
	int min1 = INT_MAX;
	vector<int> pos;
	for (int i = 0; i < s.length(); i++)
		if (set.find(s[i]) != set.end())
			pos.push_back(i);
	if (set.empty())
		return "";
	if (pos.empty())
		return "";
	while (right < pos.size()) {
		while (right < pos.size() && !LeetCode76_isdesirable(set, pos, s, left, right))
			++right;
		if (right >= pos.size())
			break;
		if (min1 > pos[right] - pos[left] + 1) {
			min1 = pos[right] - pos[left] + 1;
			minl = left;
			minr = right;
		}
		while (left + 1 <= right && LeetCode76_isdesirable(set, pos, s, left + 1, right))
			++left;
		if (min1 > pos[right] - pos[left] + 1) {
			min1 = pos[right] - pos[left] + 1;
			minl = left;
			minr = right;
		}
		++left;
	}
	if (minl != -1 && minr != -1)
		return s.substr(pos[minl], pos[minr] - pos[minl] + 1);
	else
		return "";
}
string LeetCode76_minWindow_fromweb(string S, string T) {
	if (T.size() > S.size()) return "";
	string res = "";
	int left = 0, count = 0, minLen = S.size() + 1;
	unordered_map<char, int> m;
	for (int i = 0; i < T.size(); ++i) {
		if (m.find(T[i]) != m.end()) ++m[T[i]];
		else m[T[i]] = 1;
	}
	for (int right = 0; right < S.size(); ++right) {
		if (m.find(S[right]) != m.end()) {
			--m[S[right]];
			if (m[S[right]] >= 0) ++count;//可能在S中找到2个a，在T中只有一个a
			while (count == T.size()) {
				if (right - left + 1 < minLen) {
					minLen = right - left + 1;
					res = S.substr(left, minLen);
				}
				if (m.find(S[left]) != m.end()) {
					++m[S[left]];
					if (m[S[left]] > 0) --count;//将属于T的字符取出，需要更新count
				}
				++left;
			}
		}
	}
	return res;
}
int LeetCode91_recu(string s) {
	if (s.empty())
		return 1;
	if (s[0] == '0')
		return 0;
	int num = 0;
	int i = 0;
	int count = 0;
	while(i<s.length()) {
		num = num * 10 + s[i] - '0';
		if (num > 0 && num <= 26)
			count+=LeetCode91_recu(s.substr(i+1, s.length() - i - 1));
		if (num > 26)
			break;
		i++;
	}
	return count;
}
int LeetCode91_numDecodings(string s) {
	if (s.empty())
		return 0;
	return LeetCode91_recu(s);
}
vector<int> LeetCode89_grayCode(int n)
{
	vector<int> result(1, 0);
	for (int i = 0; i < n; i++) {
		int curCount = result.size();
		// push back all element in result in reverse order
		while (curCount) {
			curCount--;
			int curNum = result[curCount];
			curNum += (1 << i);
			result.push_back(curNum);
		}
	}
	return result;
}
TreeNode* LeetCode99_first = NULL;
TreeNode* LeetCode99_second = NULL;
TreeNode* LeetCode99_prev = new TreeNode(INT_MIN);
void LeetCode99_inorder(TreeNode* root) {
	if (!root) return;
	LeetCode99_inorder(root->left);
	if (LeetCode99_first == NULL && 
		LeetCode99_prev->val > root->val) 
		LeetCode99_first = LeetCode99_prev;
	if (LeetCode99_first != NULL && 
		LeetCode99_prev->val > root->val) 
		LeetCode99_second = root;
	//在中序遍历序列中，右节点的前一个是他的父节点
	LeetCode99_prev = root;
	LeetCode99_inorder(root->right);
}
void LeetCode99_recoverTree(TreeNode* root) {
	//二叉搜索树的中序遍历是一个递增序列，被交换的两个数，一定是中序序列中，比它后面的数大的数
	//故只要在中序遍历中，比较每一个数和它的前一个即可
	LeetCode99_inorder(root);
	int tmp= LeetCode99_first->val;
	LeetCode99_first->val = LeetCode99_second->val;
	LeetCode99_second->val=tmp;
}
struct TreeLinkNode {
	int val;
	TreeLinkNode *left, *right, *next;
	TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
};
void LeetCode116_dfs(TreeLinkNode*p, int lv, vector<vector<TreeLinkNode*>>& ret) {
	if (!p)return;
	if (lv >= ret.size())
		ret.resize(ret.size() + 1);
	ret[lv].push_back(p);
	LeetCode116_dfs(p->left, lv + 1, ret);
	LeetCode116_dfs(p->right, lv + 1, ret);
}
void LeetCode116_connect(TreeLinkNode *root) {
	vector<vector<TreeLinkNode*>> ret;
	LeetCode116_dfs(root, 0, ret);
	for (int i = 0; i < ret.size(); i++) {
		if (ret[i].empty())
			continue;
		for (int j = 0; j < ret[i].size()-1; j++) 
			ret[i][j]->next = ret[i][j + 1];
		ret[i][ret[i].size() - 1]->next = NULL;
	}
}
void LeetCode116_connect_OCspace(TreeLinkNode *root) {
	TreeLinkNode dummy(0), *t = &dummy;
	while (root && root->left) {
		t = t->next = root->left;   //connect every node in the same layer starting from dummy
		t = t->next = root->right;
		if (root->next) root = root->next;
		else { root = dummy.next; t = &dummy; }  //go to the next layer
	}
}
void LeetCode117_connect(TreeLinkNode *root) {
	TreeLinkNode dummy(0), *t = &dummy;
	while (root) {
		if (root->left) t = t->next = root->left;  //connect every node in the same layer starting from dummy
		if (root->right) t = t->next = root->right;
		if (root->next) root = root->next;
		else { root = dummy.next; dummy.next = NULL; t = &dummy; }  //go to the next layer
	}
}
TreeLinkNode* LeetCode116test_MakeTree(vector<int>& arr) {
	vector<TreeLinkNode*> p;
	for (int i = 0; i < arr.size(); i++) {
		p.push_back(new TreeLinkNode(arr[i]));
		if (arr[i] != NULL && i != 0) {
			if (i % 2 == 0)
				p[(i - 1) / 2]->right = p[i];
			else
				p[(i - 1) / 2]->left = p[i];
		}
	}
	return p[0];
}
vector<vector<int>> LeetCode118_generate(int numRows) {
	vector<vector<int>> ret;
	if (numRows == 0)
		return ret;
	ret.push_back(vector<int>(1, 1));
	if (numRows == 1)
		return ret;
	ret.push_back(vector<int>(2, 1));
	if (numRows == 2)
		return ret;
	vector<int> t;
	for (int i = 2; i < numRows; i++) {
		t.clear();
		t.push_back(1);
		for (int j = 1; j <i; j++)
			t.push_back(ret[i - 1][j] + ret[i - 1][j - 1]);
		t.push_back(1);
		ret.push_back(t);
	}
	return ret;
}
vector<int> LeetCode119_getRow(int rowIndex) {
	auto t=LeetCode118_generate(rowIndex+1);
	return t[rowIndex];
}
int LeetCode120_minimumTotal(vector<vector<int>>& triangle) {
	//自底向上动归，i行j列的点只能由，i+1行j列或j+1列到达
	if (triangle.empty())
		return 0;	
	for (int i = triangle.size() - 1; i >= 0; i--) {
		for (int j = 0; j < triangle[i].size()-1; j++)
			triangle[i - 1][j] += min(triangle[i][j], triangle[i][j + 1]);
	}
	return triangle[0][0];
}
int LeetCode121_maxProfit(vector<int>& prices) {
	int minprice = INT_MAX;
	int maxprofit = 0;
	for (int i = 0; i < prices.size(); i++) {
		if (prices[i] < minprice)
			minprice = prices[i];
		else if (prices[i] - minprice > maxprofit)
			maxprofit = prices[i] - minprice;
	}
	return maxprofit;
}
int LeetCode122_maxProfit(vector<int>& prices) {
	if (prices.empty())
		return 0;
	int maxprofit = 0;
	for (int i = 0; i < prices.size()-1; i++) 
		if (prices[i + 1] > prices[i])
			maxprofit += prices[i + 1] - prices[i];
	return maxprofit;
}
bool LeetCode125_isPalindrome(string s) {
	if (s.empty())
		return true;
	int i = 0, j = s.length() - 1;
	while (i < j) {
		while (i<j && !(s[i] >= 'a' && s[i] <= 'z' || s[i] >= 'A' && s[i] <= 'Z' || s[i] >= '0' && s[i] <= '9'))
			++i;
		while (i<j && !(s[j] >= 'a' && s[j] <= 'z' || s[j] >= 'A' && s[j] <= 'Z' || s[j] >= '0' && s[j] <= '9'))
			--j;
		if(tolower(s[i]) != tolower(s[j]))
			return false;
		++i;
		--j;
	}
	return true;
}
void LeetCode129_dfs(TreeNode* root,vector<int>& nums,int num) {
	num = num * 10 + root->val;
	if (root->left)
		LeetCode129_dfs(root->left, nums,num);
	if (root->right)
		LeetCode129_dfs(root->right, nums,num);
	if (!root->left && !root->right)
		nums.push_back(num);
}
int LeetCode129_sumNumbers(TreeNode* root) {
	if (!root)
		return 0;
	vector<int> nums;
	LeetCode129_dfs(root, nums, 0);
	int count = 0;
	for (int i = 0; i < nums.size(); ++i)
		count += nums[i];
	return count;
}
int LeetCode136_singleNumber(vector<int>& nums) {
	int n = 0;
	for (int i = 0; i < nums.size(); ++i)
		n = n^nums[i];
	return n;
}
int LeetCode134_canCompleteCircuit(vector<int>& gas, vector<int>& cost) {
	int accu = 0, sind = 0,count=0,diff;
	for (int i = 0; i < gas.size(); ++i) {
		diff = gas[i] - cost[i];
		count += diff;
		if (accu + diff >= 0)
			accu += diff;
		else {
			sind = i + 1;
			accu = 0;
		}
	}
	if (count<0)
		return -1;
	else
		return sind;
}
ListNode* LeetCode148_quicksort(ListNode* h, ListNode* tn) {
	if (!h)
		return NULL;
	if (h->next == tn)
		return h;
	ListNode* p = h,
		*dum=new ListNode(-1),//设置一个辅助链表头
		*pp=dum,*ph=dum;
	dum->next = h;
	while (p!= tn) {
		if (p->val >= h->val){//以头做基准
			p = p->next;
			pp = pp->next;
		}
		else {//将小于基准的，扔到前面
			ph->next = p; 
			ph = p;
			pp->next = p->next;
			p->next = h;
			p = pp->next;
		}
	}
	pp = dum->next;
	delete dum;
	//以原来的h为分割，递归排序两边区间
	if (pp != h)//如果排序区间不为空才排序
		p = LeetCode148_quicksort(pp, h);
	else
		p = pp;
	if(h->next!=tn){
		pp=LeetCode148_quicksort(h->next, tn);
		h->next = pp;//将两个排序区间的尾头重新对接
	}
	return p;//返回此排序区间的头部
}
ListNode* LeetCode148_sortList(ListNode* head) {
	return LeetCode148_quicksort(head,NULL);
}
//将b插入到a后，bp是b原来的前一个
void LeetCode147_insert(ListNode* a, ListNode* b,ListNode* bp) {
	//ListNode* p = b->next;
	bp->next = b->next;
	b->next = a->next;
	a->next = b;
}
ListNode* LeetCode147_insertionSortList(ListNode* head) {
	if (!head)
		return NULL;
	if (!head->next)
		return head;
	int count = 1;//已排序区间长度
	ListNode* p = head->next,
		*dum=new ListNode(-1),
		*pp=head;
	dum->next = head;
	//p指向待插入节点，pp指向其前一个，a用于遍历有序区间
	while (p != NULL) {
		ListNode*a = dum->next;
		//新增节点在有序区间头插入
		if (p->val < a->val) {
			LeetCode147_insert(dum, p, pp);
			p = pp->next;
			++count;
			continue;
		}
		//在区间中间插入
		int i = 1;
		for (; i < count; ++i) {
			if (p->val >= a->val && p->val < a->next->val) {
				LeetCode147_insert(a, p, pp);
				p = pp->next;
				++count;
				break;
			}
			else 
				a = a->next;
		}
		//在区间末尾插入
		if (i == count) {
			pp = pp->next;
			p = p->next;
			++count;
		}
	}
	return dum->next;
}
struct RandomListNode {
	int label;
	RandomListNode *next, *random;
	RandomListNode(int x) : label(x), next(NULL), random(NULL) {}
};
RandomListNode *LeetCode138_copyRandomList(RandomListNode *head) {
	map<RandomListNode*, RandomListNode*> m;
	if (!head)return nullptr;
	RandomListNode*p = head,*nh,*pp;
	//新建头
	nh = new RandomListNode(head->label);
	if(head->random){
		nh->random=new RandomListNode(head->random->label);
		m[head->random] = nh->random;}
	pp = nh;
	while (p->next) {
		pp->next = new RandomListNode(p->next->label);
		p = p->next;
		pp = pp->next;
		if (!p->random)continue;
		if (m.find(p->random) != m.end()) {
			pp->random = m[p->random];
		}
		else {
			pp->random = new RandomListNode(p->random->label);
			m[p->random] = pp->random;
		}
	}
	return nh;
}
bool LeetCode139_recur(string s, vector<string>& wordDict) {
	string t1=s.substr(0,1),t2=t1;
	t2[0] = (char)t2[0] + 1;
	for (auto i = lower_bound(wordDict.begin(), wordDict.end(), t1); i != lower_bound(wordDict.begin(), wordDict.end(), t2); ++i) {
		if (s.length() >= i->length() && s.substr(0, i->length()) == *i) {
			if (s.length() == i->length())return true;
			else if (LeetCode139_recur(s.substr(i->length(), s.length() - i->length()), wordDict)) return true;
		}
	}
	return false;
}
bool LeetCode139_wordBreak(string s, vector<string>& wordDict) {
	if (s.empty())return true;
	sort(wordDict.begin(), wordDict.end());
	return LeetCode139_recur(s, wordDict);
}
bool LeetCode139_wordBreak_dp(string s, vector<string>& wordDict) {
	vector<bool> memory(s.length() + 1);
	memory[0] = true;
	for (int i = 0; i < s.length(); i++) {
		if (!memory[i]) {
			continue;
		}
		for (string word : wordDict) {
			int end = i + word.length();
			if (end <= s.length() && s.substr(i, end-i)==word) {
				memory[end] = true;
			}
		}
	}
	return memory[s.length()];
}
//LeetCode146
class LRUCache {
public:
	LRUCache(int capacity) {
		this->capacity = capacity;
	}
	int get(int key) {
		if (m.find(key) != m.end()){
			++count[key];
			history.push_back(key);
			return m[key];}
		else return -1;
	}
	void put(int key, int value) {
		if (m.find(key) != m.end()) {
			history.push_back(key);
			m[key] = value;
			++count[key];
		}
		else {//新key
			//有空余
			if (capacity > 0) {
				--capacity;
				m[key] = value;
				count[key]=1;
				history.push_back(key);
			}
			else {//无空余，按访问历史更新引用计数
				while (history.size() > 1 && count[history.front()]>1) {
					--count[history.front()];
					history.pop_front();
				}
				//删掉最少访问的
				m.erase(history.front());
				count.erase(history.front());
				history.pop_front();
				m[key] = value;
				count[key] = 1;
				history.push_back(key);
			}
		}
	}
private:
	deque<int> history;
	unordered_map<int,int> m;
	unordered_map<int,int> count;
	int capacity;
};
//LeetCode155
class MinStack {
public:
	MinStack():size(0),m(2147483647) {

	}
	void push(int x) {
		arr[size++] = x;
		if (x < m)m = x;
	}
	void pop() {
		--size;
		if (arr[size] == m) {
			if (size != 0)
				m = arr[0];
			else
				m = 2147483647;
			for (int i = 1; i < size; ++i)if (arr[i] < m)m = arr[i];
		}
	}
	int top() {
		return arr[size-1];
	}
	int getMin() {
		return m;
	}
private:
	int arr[102400];
	int size;
	int m;
};
ListNode* LeetCode160_getIntersectionNode(ListNode *headA, ListNode *headB) {
	ListNode *pa = headA, *pb = headB;
	bool a = false, b = false;
	if (!pa || !pb)return nullptr;
	while (pa != pb) {
		pa = pa->next;
		if (!pa)
			if (a) return nullptr; 
			else{
				pa = headB; a = true;}
		pb = pb->next;
		if (!pb) 
			if (b) return nullptr;
			else {
				pb = headA; b = true;}
	}
	return pa;
}
vector<int> LeetCode167_twoSum(vector<int>& numbers, int target) {
	vector<int> ret;
	for (int i = 0; i < numbers.size(); ++i) {
		int t = target - numbers[i];
		if (binary_search(numbers.begin() + i + 1, numbers.end(), t)) {
			int j = lower_bound(numbers.begin() + i + 1, numbers.end(), t) - numbers.begin();
			ret.push_back(i+1); ret.push_back(j+1);
			return ret;
		}
	}
	return ret;
}
string LeetCode168_convertToTitle(int n) {
	string ret;
	while (n != 0) {
		if (n % 26)
			{ret.insert(ret.begin(), n % 26 - 1 + 'A'); n = n / 26;}
		else
		{
			ret.insert(ret.begin(), 'Z'); n = (n - 1) / 26;}
		
	}
	return ret;
}
int LeetCode171_titleToNumber(string s) {
	int ret = 0;
	for (int i = 0; i < s.size(); ++i) 
		ret+=(s[i]-'A'+1)*pow(26, s.size() - i - 1);
	return ret;
}
int LeetCode172_trailingZeroes(int n) {
	/*要求末尾有多少个零，则该数应为x*10k 的形式等于x*（2k *5k）
	也就是求该数分解质因子后有几个5就行，
	：如1*2*3*4*5=1*2*3*2*2*5（里面有一个5）所以结果为1个0*/
	int sum = 0;
	while (n>0) {
		sum += n / 5;
		n /= 5;
	}
	return sum;
}
void LeetCode189_rotate(vector<int>& nums, int k) {
	for (int i = 0; i < k; ++i) {
		nums.insert(nums.begin(), nums.back());
		nums.pop_back();
	}
}
uint32_t LeetCode190_reverseBits(uint32_t n) {
	uint32_t mask = 0x80000000,ret=0;
	for (int i = 1; i <= 16; ++i) {
		ret |= (mask&n) >> (32 - 2*i+1);
		mask >>= 1;
	}
	for (int i = 1; i <= 16; ++i) {
		ret |= (mask&n) << (2*i-1);
		mask >>= 1;
	}
	return ret;
}
int LeetCode115_numDistinct(string s, string t) {
	vector<int> stack1;
	int count = 0;
	int i = 0;
	while (!(stack1.empty() && i >= s.length())) {
		if (i >= s.length()) {
			i = stack1.back() + 1;
			stack1.pop_back();
			continue;
		}
		if (t[stack1.size()] == s[i]){
			stack1.push_back(i);
			if(stack1.size()==t.length()){
				++count;
				i = stack1.back();
				stack1.pop_back();
			}}
		++i;
	}
	return count;
}
int LeetCode115_numDistinct_dp(string s, string t) {
	const int szS = s.size();
	const int szT = t.size();
	if (szS == 0)    return 0;
	if (szT == 0)    return 1;
	vector<long long> dp(szT + 1, 0);
	dp[0] = 1;
	for (int row = 0; row<szS; ++row){
		for (int col = szT - 1; col >= 0; --col){
			if (s[row] == t[col])  
				dp[col + 1] += dp[col];
			/*dp[0-s,0,t]
			dp[0,*]=0,dp[1,*]=1
			if(s[i]==t[j]) dp[i,j]=dp[i-1,j]+dp[i-1,j-1]
			else dp[i,j]=dp[i-1,j]
			*/
		}}
	return dp[szT];
}
int LeetCode124_maxToLeaf(TreeNode* root) {
	//if (!root)return (-2147483647 - 1)/*INT_MIN*/;
	if (!root->left && !root->right)return root->val;
	if (root->left && !root->right)return LeetCode124_maxToLeaf(root->left);
	if (!root->left && root->right)return LeetCode124_maxToLeaf(root->right);
	return max(LeetCode124_maxToLeaf(root->left), LeetCode124_maxToLeaf(root->right))+root->val;
}
int LeetCode124_maxPathSum(TreeNode* root) {
	if (!root)return (-2147483647 - 1)/*INT_MIN*/;
	if (!root->left && !root->right)return root->val;
	if (root->left && !root->right)return max(root->val,max(LeetCode124_maxPathSum(root->left),
		LeetCode124_maxToLeaf(root->left)+root->val));
	if (!root->left && root->right)return max(root->val,max(LeetCode124_maxPathSum(root->right),
		LeetCode124_maxToLeaf(root->right) + root->val));
	return max(root->val,max(LeetCode124_maxToLeaf(root->left) + LeetCode124_maxToLeaf(root->right) + root->val
		,max(LeetCode124_maxPathSum(root->left)
			,LeetCode124_maxPathSum(root->right))));
}
int LeetCode124_max = (-2147483647 - 1);
int LeetCode124_helper(TreeNode* root) {
	if (!root)
		return 0;
	int left = LeetCode124_helper(root->left);
	int right = LeetCode124_helper(root->right);
	LeetCode124_max = max(LeetCode124_max, root->val + left + right);
	int current = root->val + max(left, right);
	return (current > 0) ? current : 0;
}
int LeetCode124_maxPathSum_2(TreeNode* root) {
	LeetCode124_helper(root);
	return LeetCode124_max;
}
bool LeetCode127_oneletterdiff(const string& w1, const string& w2) {
	if (w1.length() != w2.length())return false;
	int count = 0;
	for (int i = 0; i < w1.length(); ++i)
		if (w1[i] != w2[i])++count;
	return count == 1;
}
int LeetCode127_ladderLength(string beginWord, string endWord, vector<string>& wordList) {
	if (find(wordList.begin(), wordList.end(), endWord) == wordList.end())return 0;
	wordList.push_back(beginWord);
	vector<vector<bool>> graph(wordList.size(), vector<bool>(wordList.size(), 0));
	for (int i = 0; i < graph.size(); ++i)
		for (int j = i; j < graph.size(); ++j) 
			if (i == j || LeetCode127_oneletterdiff(wordList[i], wordList[j])){
				graph[i][j] = true; graph[j][i] = true;}
	set<int> viewed;
	deque<int> queue;
	viewed.insert(wordList.size() - 1);
	queue.push_back(wordList.size() - 1);
	int count = 1;
	while (!queue.empty()) {
		int t = queue.size();
		for (int j = 0; j < t;++j){
			for (int i = 0; i < wordList.size(); ++i)
				if (graph[queue.front()][i] && viewed.find(i) == viewed.end()){
					if (wordList[i] == endWord) {
						return ++count;
					}
					queue.push_back(i);
					viewed.insert(i);}
			queue.pop_front();
		}
		++count;
	}
	return 0;
}
vector<vector<string>> LeetCode126_findLadders(string beginWord, string endWord, vector<string>& wordList) {
	vector<vector<string>> ret;
	if (find(wordList.begin(), wordList.end(), endWord) == wordList.end())return ret;
	auto t = find(wordList.begin(), wordList.end(), beginWord);
	if (t != wordList.end())wordList.erase(t);
	wordList.push_back(beginWord);
	vector<vector<bool>> graph(wordList.size(), vector<bool>(wordList.size(), 0));
	for (int i = 0; i < graph.size(); ++i)
		for (int j = i+1; j < graph.size(); ++j)
			if (LeetCode127_oneletterdiff(wordList[i], wordList[j])) {
				graph[i][j] = true; graph[j][i] = true;
			}
	//set<int> viewed;
	vector<pair<int,int>> queue;
	//viewed.insert(wordList.size() - 1);
	queue.push_back(pair<int,int>(-1,wordList.size() - 1));
	int count = 1;
	int start = 0,end=1;
	bool found = false;
	while (1) {
		if (found) {
			for (int j = start; j < end; ++j) {
				vector<string> tmp;
				if (wordList[queue[j].second] == endWord) {
					int i = j;
					while (i != -1) {
						tmp.push_back(wordList[queue[i].second]);
						i = queue[i].first;}
					reverse(tmp.begin(), tmp.end());
					ret.push_back(tmp);
				}
			}
			break;
		}
		if (/*viewed.size() >= wordList.size() ||*/ start >= end)
			break;
		for (int j = start; j < end; ++j) {
			for (int i = 0; i < wordList.size(); ++i)
				if (graph[queue[j].second][i] /*&& (viewed.find(i) == viewed.end()|| wordList[i] == endWord)*/) {
					if (wordList[i] == endWord)
						found = true;
					queue.push_back(pair<int, int>(j, i));
					//viewed.insert(i);
				}
		}
		start = end;
		end = queue.size();
		++count;
	}
	return ret;
}
void LeetCode126_findChildren(string word, unordered_set<string>& next, unordered_set<string>& dict, unordered_map<string, vector<string>>& children) {
	string parent = word;
	for (int i = 0; i < word.size(); i++) {
		char t = word[i];
		for (int j = 0; j < 26; j++) {
			word[i] = 'a' + j;
			if (dict.find(word) != dict.end()) {
				next.insert(word);
				children[parent].push_back(word);
			}
		}
		word[i] = t;
	}
}
void LeetCode126_genLadders(string beginWord, string endWord, unordered_map<string, vector<string>>& children, vector<string>& ladder, vector<vector<string>>& ladders) {
	if (beginWord == endWord) {
		ladders.push_back(ladder);
	}
	else {
		for (string child : children[beginWord]) {
			ladder.push_back(child);
			LeetCode126_genLadders(child, endWord, children, ladder, ladders);
			ladder.pop_back();
		}
	}
}
vector<vector<string>> LeetCode126_findLadders_BfsDfs(string beginWord, string endWord, vector<string>& wordList) {
	unordered_set<string> dict(wordList.begin(), wordList.end()), current, next;
	if (dict.find(endWord) == dict.end()) {
		return{};
	}
	unordered_map<string, vector<string>> children;
	vector<vector<string>> ladders;
	vector<string> ladder;
	current.insert(beginWord);
	ladder.push_back(beginWord);
	while (true) {
		for (string word : current) {
			dict.erase(word);
		}
		for (string word : current) {
			LeetCode126_findChildren(word, next, dict, children);
		}
		if (next.empty()) {
			break;
		}
		if (next.find(endWord) != next.end()) {
			LeetCode126_genLadders(beginWord, endWord, children, ladder, ladders);
			break;
		}
		current.clear();
		swap(current, next);
	}
	return ladders;
}
class Node {
public:
	int val;
	vector<Node*> neighbors;
	Node() {}
	Node(int _val, vector<Node*> _neighbors) {
		val = _val;
		neighbors = _neighbors;
	}
};
Node* LeetCode133_cloneGraph(Node* node) {
	map<Node*, Node*> m;
	if (!node)return nullptr;
	deque<Node*> que;
	que.push_back(node);
	while (!que.empty()) {
		Node*t = que.front();
		Node* nt;
		if (m.find(t) == m.end()) {
			nt = new Node();
			m[t] = nt;
			nt->val = t->val;
		}
		else
			nt = m[t];
		for (Node* p : t->neighbors) {
			Node* np;
			if (m.find(p) == m.end()) {
				np = new Node();
				m[p] = np;
				np->val = p->val;
				que.push_back(p);
			}
			else
				np = m[p];
			nt->neighbors.push_back(np);
		}
		que.pop_front();
	}
	return m[node];
}
void LeetCode130_solve(vector<vector<char>>& board) {
	if (board.empty())return;
	deque<pair<int, int>> que;
	set<pair<int, int>> set1;
	for (int i = 0;i < board.size();++i) {
		if (board[i][0] == 'O') {
			que.push_back(pair<int, int>(i, 0));
			set1.insert(pair<int, int>(i, 0));
		}
		if (board[i][board[0].size() - 1] == 'O') {
			que.push_back(pair<int, int>(i, board[0].size() - 1));
			set1.insert(pair<int, int>(i, board[0].size() - 1));
		}
	}
	for (int i = 1;i < (int)board[0].size()-1;++i) {
		if (board[0][i] == 'O') {
			que.push_back(pair<int, int>(0, i));
			set1.insert(pair<int, int>(0, i));
		}
		if (board[board.size() - 1][i] == 'O') {
			que.push_back(pair<int, int>(board.size() - 1, i));
			set1.insert(pair<int, int>(board.size() - 1, i));
		}
	}
	int m = board.size(), n = board[0].size();
	while (!que.empty()) {
		pair<int, int> t = que.front();
		if (t.first > 0 && board[t.first - 1][t.second] == 'O'&& set1.find(pair<int, int>(t.first - 1, t.second)) == set1.end()) {
			set1.insert(pair<int, int>(t.first - 1, t.second));
			que.push_back(pair<int, int>(t.first - 1, t.second));
		}
		if (t.first < m-1 && board[t.first + 1][t.second] == 'O'&& set1.find(pair<int, int>(t.first + 1, t.second)) == set1.end()) {
			set1.insert(pair<int, int>(t.first + 1, t.second));
			que.push_back(pair<int, int>(t.first + 1, t.second));
		}
		if (t.second > 0 && board[t.first ][t.second - 1] == 'O'&& set1.find(pair<int, int>(t.first, t.second - 1)) == set1.end()) {
			set1.insert(pair<int, int>(t.first , t.second - 1));
			que.push_back(pair<int, int>(t.first , t.second - 1));
		}
		if (t.second < n-1 && board[t.first][t.second+1] == 'O'&& set1.find(pair<int, int>(t.first, t.second + 1))==set1.end()) {
			set1.insert(pair<int, int>(t.first , t.second+1));
			que.push_back(pair<int, int>(t.first , t.second+1));
		}
		que.pop_front();
	}
	for (int i = 0;i < m;++i)
		for (int j = 0;j < n;++j) {
			if (board[i][j] == 'O'&& set1.find(pair<int, int>(i, j) )== set1.end())
				board[i][j] = 'X';
		}
}
int LeetCode128_longestConsecutive(vector<int>& nums) {
	unordered_set<int> set1;
	for (int num : nums)set1.insert(num);
	int len = 0;
	for (int num : nums) {
		int l = 1;
		int cn = num;
		while (set1.find(cn + 1) != set1.end()) {
			++l;
			++cn;
		}
		len = max(l, len);
	}
	return len;
}
bool LeetCode131_isPalindrome(const string &s, int i, int j) {
	while (i < j) {
		if (s[i] != s[j])return false;
		++i;
		--j;
	}
	return true;
}
int LeetCode131_findNextAvailable(const string &s, const pair<int, int>& pa) {
	int j = pa.second + 1;
	for (;j < s.size();++j)
		if (LeetCode131_isPalindrome(s, pa.first, j))return j;
	return -1;
}
vector<vector<string>> LeetCode131_partition(string s) {
	vector<pair<int, int>> stack1;
	vector<vector<string>> ret;
	int l;
	stack1.push_back(pair<int, int>(0, 0));
	while (!stack1.empty()){
		//到达串尾，输出结果
		if (stack1.back().second == s.size() - 1) {
			vector<string> tmp;
			for (auto t : stack1) 
				tmp.push_back(s.substr(t.first, t.second - t.first + 1));
			ret.push_back(tmp);
			//弹出栈并搜索下一个可用分割
			while(1){
				stack1.pop_back();
				if (stack1.empty())return ret;
				int j = LeetCode131_findNextAvailable(s,stack1.back());
				if (j != -1) {
					stack1.back().second = j;
					break;
				}
			}
		}
		pair<int, int> t(stack1.back().second + 1, stack1.back().second);
		int j = LeetCode131_findNextAvailable(s, t);
		if (j != -1) {
			t.second = j;
			stack1.push_back(t);
		}
	}
	return ret;
}
int LeetCode132_minCut(string s) {
	if (s.empty())return 0;
	int current = 1;
	deque<int> curque,nextque;
	curque.push_back(-1);
	while (!curque.empty()) {
		int t = curque.front() + 1;
		for (int i = t;i < s.size();++i)
			if (LeetCode131_isPalindrome(s, t, i)) {
				if (i == s.size() - 1)return current - 1;
				nextque.push_back(i);
			}
		curque.pop_front();
		if (curque.empty()) {
			swap(curque, nextque);
			++current;
		}
	}
	return -1;
}
int LeetCode132_minCut_dp(string s) {
	int n = s.length();
	vector<vector<bool>> dp(n,vector<bool>(n,false));
	vector<int> cut(n,n-1);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i <= j; i++) {
			if (s[i] == s[j] && (j - i <= 1 || dp[i + 1][j - 1])) {
				dp[i][j] = true;
				if (i == 0)
					cut[j] = 0; //no cut needed
				else
					cut[j] = min(cut[j], cut[i - 1] + 1);
			}
		}
	}
	return cut[n - 1];
}
int LeetCode137_singleNumber(vector<int>& nums) {
	int bitSum[32] = { 0 };
	int n = nums.size();
	int result = 0;
	for (int i = 0; i < 32; i++) {
		for (int j = 0; j < n; j++) {
			bitSum[i] += (nums[j] >> i) & 1;
		}
		bitSum[i] %= 3;
		result |= bitSum[i] << i;
	}
	return result;
}
int main() {
	auto t =LeetCode132_minCut_dp("abbab");
	int b = 1;
}
