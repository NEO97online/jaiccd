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

CCD_List :: struct {
    next, prev: *CCD_List;
}

ccd_offsetof :: ($T: Type, ident: Code) -> s64 #expand {
    t: T = ---;
    return cast(*void) (*t.#insert ident) - cast(*void) *t;
}

/**
 * Get the struct for this entry.
 * @ptr:	the &CCD_List pointer.
 * @type:	the type of the struct this is embedded in.
 * @member:	the name of the list_struct within the struct.
 */
ccdListEntry :: (ptr: *CCD_List, $T: Type) -> *T {
    return cast(*T) ( cast(*void)(ptr) - ccd_offsetof(T, #code list) );
}

/**
 * Iterates over list of given type.
 * @pos:	the type * to use as a loop cursor.
 * @head:	the head for your list.
 * @member:	the name of the list_struct within the struct.
 */
ccdListForEachEntry :: (head: *CCD_List, pos: **PosType, $PosType: Type, code: Code) #expand {
    pos.* = ccdListEntry(head.next, PosType);
    while *(pos.*.list) != head {
        #insert code;
        pos.* = ccdListEntry(pos.*.list.next, PosType);
        _ccd_prefetch(pos.*.list.next);
    }
}

///
/// INLINES:
///

/**
 * Initialize list.
 */
ccdListInit :: inline (l: *CCD_List) {
    l.next = l;
    l.prev = l;
}

ccdListNext :: inline (l: *CCD_List) -> *CCD_List {
    return l.next;
}

ccdListPrev :: inline (l: *CCD_List) -> *CCD_List {
    return l.prev;
}

/**
 * Returns true if list is empty.
 */
ccdListEmpty :: inline (head: *CCD_List) -> bool {
    return head.next == head;
}

/**
 * Appends item to end of the list l.
 */
ccdListAppend :: inline (l: *CCD_List, new: *CCD_List) {
    new.prev = l.prev;
    new.next = l;
    l.prev.next = new;
    l.prev = new;
}

/**
 * Removes item from list.
 */
ccdListDel :: inline (item: *CCD_List) {
    item.next.prev = item.prev;
    item.prev.next = item.next;
    item.next = item;
    item.prev = item;
}
