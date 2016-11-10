#include <iostream>
#include <vector>

struct ListNode {
  int val;
  ListNode *next;
  ListNode(int x) : val(x), next(NULL) {}
};
ListNode* reverseBetween(ListNode* head, int m, int n);
void nextStep(ListNode*&, ListNode*&, ListNode*&);
void listSwap(ListNode*&, ListNode*&, ListNode*&);

int main() {
  ListNode *head = NULL, *temp = NULL, *iter = NULL;
  std::vector<int> vector;

  //Create the list
  for (int i=0; i<10; ++i) {
    vector.push_back(i);
  }
  for (int i=vector.size()-1; i>=0; --i) {
    ListNode* newNode = new ListNode(vector[i]);
    newNode->next = temp;
    temp = newNode;
  }
  head = temp;

  head = reverseBetween(head, 1, 1);

  iter = head;
  while (iter) {
    std::cout << iter->val << " ";
    iter = iter->next;
  }
  std::cout << std::endl;

  //Clean up the list
  while (temp) {
    ListNode* next = temp->next;
    delete temp;
    temp = next;
  }
}

/* reverses*/
void listSwap(ListNode*& prev, ListNode*& curr, ListNode*& forward) {
  curr->next = prev;
  nextStep(prev, curr, forward);
}

/* Re-assigns PREV, CURR, and FORWARD to be the next sequence of three in the linked list*/
void nextStep(ListNode*& prev, ListNode*& curr, ListNode*& forward) {
  prev = curr;
  curr = forward;
  if (forward) forward = forward->next;
}

ListNode* reverseBetween(ListNode* head, int m, int n) {
  if (m==n || !head->next) return head;    //no reversing necessary

  ListNode* newHead = head;
  ListNode* start = NULL; //element at the start of the flippy section
  ListNode* leftBound = NULL; //element immediately before the flippy section

  ListNode* backward = NULL;
  ListNode* current = head;
  ListNode* forward = head->next;
  int i=1;

  while (current && i<=n) {
    if (i==m) {leftBound = backward; start = current;}

    if (i>m) {
      listSwap(backward, current, forward);
    }
    else {
      nextStep(backward, current, forward);
    }
    ++i;
  }

  if (!leftBound) newHead = backward;
  else leftBound->next = backward;
  start->next = current;

  return newHead;
 }
