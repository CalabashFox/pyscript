import re

references_dict = dict()
misses = []


class Data(object):
    def __init__(self, index, gen):
        self.index = index
        self.gen = gen
        self.gen_seq = 0
        self.gen_id = 0
        self.fusions = []
        self.tissues = []
        self.frameshifts = dict()
        self.snvs = dict()


class Fusion(object):
    def __init__(self, mutation, mutation_id, count, link):
        self.mutation = mutation
        self.mutation_id = mutation_id
        self.count = count
        self.link = link

    def __repr__(self):
        return 'mutation: {}, mutation_id: {}, count: {}, link: {}'\
            .format(self.mutation, self.mutation_id, self.count, self.link)


class Tissue(object):
    def __init__(self, variation_type, conseq, value, total, positive):
        self.variation_type = variation_type
        self.conseq = conseq
        self.value = value
        self.total = total
        self.positive = positive

    def __repr__(self):
        return 'type: {}, conseq: {}, value: {}, total: {}, positive: {}'\
            .format(self.variation_type, self.conseq, self.value, self.total, self.positive)


mut_pattern = re.compile('\d+')


class Frameshift(object):
    def __init__(self, mut_id, transcript, variation_type, conseq, aa_mut, cds_mut, link, references):
        self.type = 0
        self.mut_id = mut_id
        self.transcript = transcript
        self.variation_type = variation_type
        self.conseq = conseq
        self.aa_mut = aa_mut
        self.cds_mut = cds_mut
        self.link = link
        self.count = 1
        self.references = references
        n = re.findall('\d+', self.aa_mut)
        if len(n) == 0:
            self.order_mut = 0
        else:
            self.order_mut = int(n[0])

    def incr(self):
        self.count += 1

    def __repr__(self):
        return 'transcript: {}, type: {}, conseq: {}, aa_mut: {}, cds_mut: {}, link: {}'\
            .format(self.transcript, self.variation_type, self.conseq, self.aa_mut, self.cds_mut, self.link)


class SNV(object):
    def __init__(self, mut_id, transcript, conseq, aa_mut, cds_mut, score, link, references):
        self.type = 1
        self.mut_id = mut_id
        self.transcript = transcript
        self.conseq = conseq
        self.aa_mut = aa_mut
        self.cds_mut = cds_mut
        self.score = score
        self.link = link
        self.count = 1
        self.references = references
        n = re.findall('\d+', self.aa_mut)
        if len(n) == 0:
            self.order_mut = 0
        else:
            self.order_mut = int(n[0])

    def incr(self):
        self.count += 1

    def __repr__(self):
        return 'transcript: {}, conseq: {}, aa_mut: {}, cds_mut: {}, score: {}, link: {}'\
            .format(self.transcript, self.conseq, self.aa_mut, self.cds_mut, self.score, self.link)
