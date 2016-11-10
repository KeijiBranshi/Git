/********************************************************************************************
 * Palindrome detector for a linked list of ints. Accomplished in O(n) time with O(1) space.
 * Description of my method is written above the isPalindrome() function.
 ********************************************************************************************/

#include <iostream>
#include <vector>

struct ListNode {
  int val;
  ListNode* next;
  ListNode(int x) : val(x), next(NULL) {}
};

int getListSize(ListNode* head);
ListNode* getTail(ListNode* head);
ListNode* reverseList(ListNode* head);
bool isPalindrome(ListNode* head);

int main() {
  ListNode *head = NULL, *temp = NULL;
  std::vector<int> vector(1,1);

  //Create the list
  for (int i=vector.size()-1; i>=0; --i) {
    ListNode* newNode = new ListNode(vector[i]);
    newNode->next = temp;
    temp = newNode;
  }
  head = temp;

  std::cout << std::boolalpha << isPalindrome(head) << std::endl;

  //Clean up the list
  while (temp) {
    ListNode* next = temp->next;
    delete temp;
    temp = next;
  }
}

/*
 * Returns true if a linked list of ints, starting at HEAD, is a palindrome. Runs in O(n) time with O(1) space complexity.
 * My method was to split the linked list in half, reverse the second half, then make side by side comparisons
 * of the two new linked lists. This only incurs 4 O(n) loops through the list, plus an additional 2 O(n) loops
 * to restore the linked list back to normal. The last 2 loops are for cleanup purposes, since my list is
 * dynamically allocated.
 */
bool isPalindrome(ListNode* head) {
  if (!head || !head->next) return true;

  ListNode* listFirstHead = head;
  ListNode* listFirstIter = head;
  ListNode* listSecondHead = NULL;
  ListNode* listSecondIter = NULL;
  int size = getListSize(head);

  //SPLIT THE LIST
  for (int i=0; i<size/2-1; ++i){
    listFirstIter = listFirstIter->next;
  }
  listSecondHead = listFirstIter->next; //found second half
  listSecondHead = reverseList(listSecondHead);
  listFirstIter->next = NULL;   //cuts first half off from second half
  //Note that second half of linked list will alway be 0 or 1 element longer than the first half

  //ITERATE THROUGH THE TWO HALFS TO CHECK IF THE LINKED LIST IS A PALINDROME
  listFirstIter = listFirstHead;
  listSecondIter = listSecondHead;
  while (listFirstIter && listSecondIter) {
    if (listFirstIter->val != listSecondIter->val) return false;
    else {
      listFirstIter = listFirstIter->next;
      listSecondIter = listSecondIter->next;
    }
  }

  //Optional-ish: Restore the Linked List (for cleanup later)
  listSecondHead = reverseList(listSecondHead);
  listFirstIter = getTail(listFirstHead);
  listFirstIter->next = listSecondHead;

  return true;
}

//Returns the size of the linked list starting at HEAD
int getListSize(ListNode* head) {
  ListNode* curr = head;
  int size = 0;
  while (curr) {
    ++size; curr = curr->next;
  }
  return size;
}

// Returns the tail of a linked list starting at HEAD
ListNode* getTail(ListNode* head) {
  if (!head) return NULL;

  ListNode* temp = head;
  while (temp->next) {
    temp = temp->next;
  }
  return temp;
}

//Reverses order of list starting at HEAD, then returns the head of the new list
ListNode* reverseList(ListNode* head) {
  ListNode* prev = NULL;
  ListNode* curr = head;
  ListNode* next = head->next;
  while (next) {
    curr->next = prev;
    prev = curr;
    curr = next;
    next = curr->next;
  }
  curr->next = prev;
  return curr;
}
