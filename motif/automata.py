from collections import deque

class NFA:
    def __init__(self, num_states, start, accepts):
        self.num_states = num_states
        self.start = start
        self.accepts = set(accepts)
        self.delta = {}  # (state, symbol) -> set(next_states)
    def add_transition(self, s, symbol, t):
        self.delta.setdefault((s, symbol), set()).add(t)

def build_linear_nfa(symbol_sets):
    m = len(symbol_sets)
    nfa = NFA(num_states=m+1, start=0, accepts={m})
    for i, allowed in enumerate(symbol_sets):
        for sym in allowed:
            nfa.add_transition(i, sym, i+1)
    return nfa

class DFA:
    def __init__(self, start, accepts):
        self.start = start
        self.accepts = set(accepts)
        self.delta = {}  # (frozenset_state, symbol) -> frozenset_state

def nfa_to_dfa(nfa, alphabet):
    start = frozenset({nfa.start})
    accepts = set()
    if any(s in nfa.accepts for s in start):
        accepts.add(start)
    dfa = DFA(start, accepts)
    q = deque([start])
    seen = {start}
    while q:
        S = q.popleft()
        for a in alphabet:
            move = set()
            for s in S:
                move |= nfa.delta.get((s, a), set())
            T = frozenset(move)
            if not T:
                continue
            dfa.delta[(S, a)] = T
            if T not in seen:
                seen.add(T)
                if any(t in nfa.accepts for t in T):
                    dfa.accepts.add(T)
                q.append(T)
    return dfa
