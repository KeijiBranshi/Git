#ifndef EVENT_H
#define EVENT_H

class Event {
protected:
  double rate_;
public:
  static const NullEvent NULLEVENT;
  double delay_;

  Event(double rate);

  double getRate();
  virtual void processEvent(TransitionMatrix& tm) = 0;
  virtual void renderEvent() = 0;
};

class InjectionEvent : public Event {
private:
  int hemeIdx_;
  Heme* heme;

public:
  InjectionEvent(double rate, int heme);

  virtual void processEvent(TransitionMatrix& tm);
  virtual void renderEvent();
};

class EjectionEvent : public Event {
private:
  int hemeIdx;
  Heme* heme;

public:
  EjectionEvent(double rate, int heme);

  virtual void processEvent(TransitionMatrix& tm);
  virtual void renderEvent();
};

class TransferEvent : public Event {
private:
  int from_hemeIdx;
  int to_hemeIdx;
  Heme* from_heme;
  Heme* to_heme;

public:
  TransferEvent(double rate, int from_heme, int to_heme);

  virtual void processEvent(TransitionMatrix& tm);
  virtual void renderEvent();
};

class NullEvent : public Event {
public:
  NullEvent();

  virtual void processEvent(TransitionMatrix& tm);
  virtual void renderEvent();
};

#endif
