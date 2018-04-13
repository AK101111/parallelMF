#ifndef SCHED_H_
#define SCHED_H_

typedef struct block{
  int x_index;
  int y_index;
  matrix_size MB;
  int free; // change to mutex??
  int num_updates;
  int priority;
}block;

typedef priorityQ std::priority_queue<block, std::vector<block>, lessThan>

struct lessThan{
  bool operator()(const block& lhs, const block& rhs) const{
    return lhs.priority < rhs.priority;
  }
};

void init_sched(priorityQ& pq);
block get_block();
void push_block(block b, priorityQ& pq);

#endif // SCHED_H_