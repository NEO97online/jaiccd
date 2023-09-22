/***
 * libccd
 * ---------------------------------
 * Copyright (c)2010 Daniel Fiser <danfis@danfis.cz>
 *
 *
 *  This file is part of libccd.
 *
 *  Distributed under the OSI-approved BSD License (the "License");
 *  see accompanying file BDS-LICENSE for details or see
 *  <http://www.opensource.org/licenses/bsd-license.php>.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the License for more information.
 */

ccd_list_t :: struct {
    next, prev: *ccd_list_t;
}

/**
 * Get the struct for this entry.
 * @ptr:	the &ccd_list_t pointer.
 * @type:	the type of the struct this is embedded in.
 * @member:	the name of the list_struct within the struct.
 */
ccdListEntry :: (ptr: *ccd_list_t, type: Type, member: $A) {
    return ccd_container_of(ptr, type, member);
}

/**
 * Iterates over list.
 */
ccdListForEach :: (list: $A, item: $B, code: Code) #expand {
    item = list.next;
    while item != list {
        #insert code;
        item = item.next;
        _ccd_prefetch(item.next);
    }
}

/**
 * Iterates over list safe against remove of list entry
 */
ccdListForEachSafe :: (list: $A, item: $B, tmp: $C, code: Code) #expand {
    item = list.next;
    tmp = item.next;
    while item != list {
        #insert code;
        item = tmp;
        tmp = item.next;
    }
}

/**
 * Iterates over list of given type.
 * @pos:	the type * to use as a loop cursor.
 * @head:	the head for your list.
 * @member:	the name of the list_struct within the struct.
 */
ccdListForEachEntry :: (head: *ccd_list_t, pos: $A, postype: $B, member: $C, code: Code) #expand {
    pos = ccdListEntry((head).next, postype, member);
    while *pos.member != head {
        #insert code;
        pos = ccdListEntry(pos.member.next, postype, member);
        _ccd_prefetch(pos.member.next);
    }
}

/**
 * Iterates over list of given type safe against removal of list entry
 * @pos:	the type * to use as a loop cursor.
 * @n:		another type * to use as temporary storage
 * @head:	the head for your list.
 * @member:	the name of the list_struct within the struct.
 */
// ccdListForEachEntrySafe :: (head, pos, postype, n, ntype, member) #expand {
//     for (pos = ccdListEntry((head).next, postype, member),
// 		 n = ccdListEntry(pos.member.next, postype, member);
// 	     *pos.member != (head);
// 	     pos = n, n = ccdListEntry(n.member.next, ntype, member))
// }


///
/// INLINES:
///

/**
 * Initialize list.
 */
ccdListInit :: inline (l: *ccd_list_t) {
    l.next = l;
    l.prev = l;
}

ccdListNext :: inline (l: *ccd_list_t) -> *ccd_list_t {
    return l.next;
}

ccdListPrev :: inline (l: *ccd_list_t) -> *ccd_list_t {
    return l.prev;
}

/**
 * Returns true if list is empty.
 */
ccdListEmpty :: inline (head: *ccd_list_t) -> int {
    return head.next == head;
}

/**
 * Appends item to end of the list l.
 */
ccdListAppend :: inline (l: *ccd_list_t, new: *ccd_list_t) {
    new.prev = l.prev;
    new.next = l;
    l.prev.next = new;
    l.prev = new;
}

/**
 * Removes item from list.
 */
ccdListDel :: inline (item: *ccd_list_t) {
    item.next.prev = item.prev;
    item.prev.next = item.next;
    item.next = item;
    item.prev = item;
}
